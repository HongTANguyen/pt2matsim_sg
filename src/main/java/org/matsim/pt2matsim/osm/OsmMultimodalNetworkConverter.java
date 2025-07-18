/*
 * *********************************************************************** *
 * project: org.matsim.*                                                   *
 *                                                                         *
 * *********************************************************************** *
 *                                                                         *
 * copyright       : (C) 2014 by the members listed in the COPYING,        *
 *                   LICENSE and WARRANTY file.                            *
 * email           : info at matsim dot org                                *
 *                                                                         *
 * *********************************************************************** *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *   See also COPYING, LICENSE and WARRANTY file                           *
 *                                                                         *
 * *********************************************************************** *
 */

package org.matsim.pt2matsim.osm;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Objects;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.matsim.api.core.v01.Id;
import org.matsim.api.core.v01.IdMap;
import org.matsim.api.core.v01.TransportMode;
import org.matsim.api.core.v01.network.Link;
import org.matsim.api.core.v01.network.Network;
import org.matsim.api.core.v01.network.Node;
import org.matsim.core.config.ConfigGroup;
import org.matsim.core.network.NetworkUtils;
import org.matsim.core.network.turnRestrictions.DisallowedNextLinks;
import org.matsim.core.network.turnRestrictions.DisallowedNextLinksUtils;
import org.matsim.core.network.turnRestrictions.TurnRestrictionsNetworkCleaner;
import org.matsim.core.utils.collections.CollectionUtils;
import org.matsim.core.utils.geometry.CoordUtils;
import org.matsim.core.utils.geometry.CoordinateTransformation;
import org.matsim.core.utils.geometry.transformations.IdentityTransformation;
import org.matsim.core.utils.geometry.transformations.TransformationFactory;
import org.matsim.pt2matsim.config.OsmConverterConfigGroup;
import org.matsim.pt2matsim.osm.LinkGeometryExporter.LinkDefinition;
import org.matsim.pt2matsim.osm.lib.AllowedTagsFilter;
import org.matsim.pt2matsim.osm.lib.Osm;
import org.matsim.pt2matsim.osm.lib.OsmData;
import org.matsim.pt2matsim.osm.lib.OsmXmlFileReader;
import org.matsim.pt2matsim.osm.lib.OsmDataImpl;
import org.matsim.pt2matsim.tools.NetworkTools;

// Import PBF related libraries
import com.slimjars.dist.gnu.trove.iterator.TLongObjectIterator;
import de.topobyte.osm4j.core.access.OsmIterator;
import de.topobyte.osm4j.core.dataset.InMemoryMapDataSet;
import de.topobyte.osm4j.core.dataset.MapDataSetLoader;
import de.topobyte.osm4j.core.model.iface.*;
import de.topobyte.osm4j.pbf.seq.PbfIterator;

import com.google.common.base.Verify;

/**
 * Converts {@link OsmData} to a MATSim network, uses a config file
 * ({@link OsmConverterConfigGroup}) to store conversion parameters and default
 * values.
 * <p>
 * See OSM wiki for more documentation on the consumed data:
 * <dl>
 * <dt>lanes</dt>
 * <dd>https://wiki.openstreetmap.org/wiki/Key:lanes</dd>
 * <dt>freespeed / maxspeed</dt>
 * <dd>https://wiki.openstreetmap.org/wiki/Key:maxspeed</dd>
 * </dl>
 *
 * @author polettif
 * @author mstraub - Austrian Institute of Technology
 */
public class OsmMultimodalNetworkConverter {

	private static final Logger log = LogManager.getLogger(OsmMultimodalNetworkConverter.class);

	private static final String OSM_TURN_RESTRICTION_ATTRIBUTE_NAME = OsmTurnRestriction.class.getSimpleName();
	private static final String OSM_SPECIAL_LANE = "_spec";

	public Network getNetwork() {
		return network;
	}

	/**
	 * mode == null means "all modes"
	 */
	record OsmTurnRestriction(@Nullable Set<String> modes, List<Id<Osm.Way>> nextWayIds,
							  @Nullable Id<Osm.Node> viaNodeId, RestrictionType restrictionType) {

		enum RestrictionType {
			PROHIBITIVE, // no_*
			MANDATORY; // only_*
		}

	}

	private static final Map<String, String> OSM_2_MATSIM_MODE_MAP = Map.of(
			Osm.Key.BUS, "bus",
			Osm.Value.TROLLEYBUS, "bus",
			Osm.Key.BICYCLE, TransportMode.bike,
			Osm.Key.MOTORCYCLE, TransportMode.motorcycle,
			Osm.Key.MOTORCAR, TransportMode.car);

	private static final List<String> TURN_RESTRICTION_KEY_SUFFIXES = List.of(
			"", // for all modes
			":" + Osm.Key.BUS,
			":" + Osm.Key.BICYCLE,
			":" + Osm.Key.MOTORCAR);

	static final int SPEED_LIMIT_WALK_KPH = 10;
	// // no speed limit (Germany) .. assume 200kph
	static final int SPEED_LIMIT_NONE_KPH = 200;

	protected final OsmData osmData;
	protected final Map<String, Map<String, OsmConverterConfigGroup.OsmWayParams>> wayParams = new HashMap<>();
	/**
	 * Maps for unknown entities
	 */
	protected final Set<String> unknownHighways = new HashSet<>();
	protected final Set<String> unknownRailways = new HashSet<>();
	protected final Set<String> unknownWays = new HashSet<>();
	protected final Set<String> unknownMaxspeedTags = new HashSet<>();
	protected final Set<String> unknownLanesTags = new HashSet<>();
	/**
	 * connects osm way ids and link ids of the generated network
	 **/
	protected final Map<Id<Link>, Id<Osm.Way>> osmIds = new HashMap<>();
	/**
	 * From one OSM way, multiple MATSim links can be created:
	 * 1) forward & reverse links
	 * 2) long or curvy ways will be split into multiple links
	 */
	protected final Map<Id<Osm.Way>, List<Id<Link>>> wayLinkMap = new HashMap<>(); // reverse of osmIds
	protected final Map<Id<Link>, DisallowedNextLinks> disallowedNextLinks = new IdMap<>(Link.class);
	protected OsmConverterConfigGroup config;
	protected Network network;
	protected long id = 0;

	protected AllowedTagsFilter ptFilter;
	protected OsmConverterConfigGroup.OsmWayParams ptDefaultParams;
	protected LinkGeometryExporter geometryExporter;

	public OsmMultimodalNetworkConverter(OsmData osmData) {
		this.osmData = osmData;
	}

	/**
	 * Converts the OSM data according to the parameters defined in config.
	 */
	public void convert(OsmConverterConfigGroup config) {
		this.config = config;
		this.geometryExporter = new LinkGeometryExporter();
		CoordinateTransformation transformation = (config.getOutputCoordinateSystem() == null ?
				new IdentityTransformation() :
				TransformationFactory.getCoordinateTransformation(TransformationFactory.WGS84, config.getOutputCoordinateSystem()));

		initPT();
		readWayParams();
		convertToNetwork(transformation);
		if (config.parseTurnRestrictions) {
			addDisallowedNextLinksAttributes();
		}
		cleanNetwork();
		if (config.getKeepTagsAsAttributes()) addAttributes();

		if (this.config.getOutputDetailedLinkGeometryFile() != null) {
			try {
				geometryExporter.onlyKeepGeometryForTheseLinks(network.getLinks().keySet());
				geometryExporter.writeToFile(Paths.get(this.config.getOutputDetailedLinkGeometryFile()));
			} catch (IOException e) {
				log.warn("Error while writing network geometry", e);
				e.printStackTrace();
			}
		}
	}

	/**
	 * reads the params from the config to different containers.
	 */
	private void readWayParams() {
		for (ConfigGroup e : config.getParameterSets(OsmConverterConfigGroup.OsmWayParams.SET_NAME)) {
			OsmConverterConfigGroup.OsmWayParams w = (OsmConverterConfigGroup.OsmWayParams) e;
			wayParams.putIfAbsent(w.getOsmKey(), new HashMap<>());
			wayParams.get(w.getOsmKey()).put(w.getOsmValue(), w);
		}
	}

	/**
	 * Converts the parsed OSM data to MATSim nodes and links.
	 */
	protected void convertToNetwork(CoordinateTransformation transformation) {

		log.info("Converting OSM to MATSim network...");

		if (transformation == null) {
			transformation = TransformationFactory.getCoordinateTransformation("WGS84", "WGS84");
		}

		this.network = NetworkTools.createNetwork();

		Map<Id<Osm.Node>, Osm.Node> nodes = osmData.getNodes();
		Map<Id<Osm.Way>, Osm.Way> ways = osmData.getWays();
		Map<Id<Osm.Relation>, Osm.Relation> relations = osmData.getRelations();

		AllowedTagsFilter serviceRailTracksFilter = new AllowedTagsFilter();
		serviceRailTracksFilter.add(Osm.ElementType.WAY, Osm.Key.SERVICE, null);

		for (Osm.Node node : nodes.values()) {
			node.setCoord(transformation.transform(node.getCoord()));
		}

		// remove ways without default params
		log.info("remove unusable ways...");
		for (Osm.Way way : new HashSet<>(ways.values())) {
			if (getWayDefaultParams(way) == null) {
				osmData.removeWay(way.getId());
			}
		}

		// remove unused nodes
		log.info("remove nodes without ways...");
		for (Osm.Node n : new HashSet<>(nodes.values())) {
			if (n.getWays().isEmpty()) {
				osmData.removeNode(n.getId());
			}
		}

		HashSet<Osm.Node> nodesToIgnore = new HashSet<>();

		log.info("cleaning network...");

		// Clean network:
		if (!config.getKeepPaths()) {
			// marked nodes as unused where only one way leads through
			// but only if this doesn't lead to links longer than MAX_LINKLENGTH
			for (Osm.Way way : ways.values()) {

				double length = 0.0;
				Osm.Node lastNode = way.getNodes().getFirst();
				for (int i = 1; i < way.getNodes().size() - 1; i++) {
					Osm.Node node = way.getNodes().get(i);
					if (node.getWays().size() > 1) {
						length = 0.0;
						lastNode = node;
					} else if (node.getWays().size() == 1) {
						length += CoordUtils.calcEuclideanDistance(lastNode.getCoord(), node.getCoord());
						if (length <= config.getMaxLinkLength()) {
							nodesToIgnore.add(node);
							lastNode = node;
						} else {
							length = 0.0;
							lastNode = node;
						}
					} else {
						log.warn("Way node with less than 1 way found.");
					}
				}
				// fix for some roundabouts with identical first and last node
				if (way.getNodes().getFirst().equals(way.getNodes().getLast())) {
					nodesToIgnore.remove(way.getNodes().getFirst());
				}
			}
			// verify we did not mark nodes as unused that build a loop
			for (Osm.Way way : ways.values()) {
				int prevRealNodeIndex = 0;
				Osm.Node prevRealNode = way.getNodes().get(prevRealNodeIndex);

				for (int i = 1; i < way.getNodes().size(); i++) {
					Osm.Node node = way.getNodes().get(i);
					if (nodesToIgnore.contains(node)) {
						if (prevRealNode == node) {
							/* We detected a loop between two "real" nodes.
							 * Set some nodes between the start/end-loop-node to "used" again.
							 * But don't set all of them to "used", as we still want to do some network-thinning.
							 * I decided to use sqrt(.)-many nodes in between...
							 */
							double increment = Math.sqrt(i - prevRealNodeIndex);
							double nextNodeToKeep = prevRealNodeIndex + increment;
							for (double j = nextNodeToKeep; j < i; j += increment) {
								int index = (int) Math.floor(j);
								Osm.Node intermediaryNode = way.getNodes().get(index);
								nodesToIgnore.remove(intermediaryNode);
							}
						}
						prevRealNodeIndex = i;
						prevRealNode = node;
					}
				}
			}
		}

		// create the required nodes and add them to the network
		log.info("Creating nodes...");
		for (Osm.Node node : nodes.values()) {
			if (!nodesToIgnore.contains(node)) {
				Node nn = this.network.getFactory().createNode(Id.create(node.getId(), Node.class), node.getCoord());
				this.network.addNode(nn);
			}
		}

		// create the links
		log.info("Creating links...");
		this.id = 1;
		int linkCount = 0;
		for (Osm.Way way : ways.values()) {
			// skip ways without nodes
			if (way.getNodes().isEmpty()) {
				log.warn("Skipping way with no nodes: {}", way.getId());
				linkCount++;
				continue;
			}
			Osm.Node fromNode = way.getNodes().getFirst();
			double length = 0.0;
			Osm.Node lastToNode = fromNode;
			if (!nodesToIgnore.contains(fromNode)) {
				for (int i = 1, n = way.getNodes().size(); i < n; i++) {
					Osm.Node toNode = way.getNodes().get(i);
					if (toNode != lastToNode) {
						length += CoordUtils.calcEuclideanDistance(lastToNode.getCoord(), toNode.getCoord());
						if (!nodesToIgnore.contains(toNode)) {
							createLink(way, fromNode, toNode, length);
							fromNode = toNode;
							length = 0.0;
						}
						lastToNode = toNode;
					}
				}
			}
		}
		log.info("Number of links created: {} and number of links having no nodes to ignore: {}",
				this.network.getLinks().size(), linkCount);

		// create reverse lookup map for link ids
		wayLinkMap.putAll(osmIds.entrySet().stream().collect(
				Collectors.groupingBy(Entry::getValue, Collectors.mapping(Entry::getKey, Collectors.toList()))));

		// parse turn restriction relations into disallowed links
		this.attachTurnRestrictionsAsDisallowedNextLinks();

		log.info("= conversion statistics: ==========================");
		log.info("MATSim: # nodes created: {}", this.network.getNodes().size());
		log.info("MATSim: # links created: {}", this.network.getLinks().size());

		if (!this.disallowedNextLinks.isEmpty()) {
			log.info("MATSim: # DisallowedNextLinks attributes created: {}", this.disallowedNextLinks.size());
		}

		if (!this.unknownHighways.isEmpty()) {
			log.info("The following highway-types had no defaults set and were thus NOT converted:");
			for (String highwayType : this.unknownHighways) {
				log.info("- \"{}\"", highwayType);
			}
		}
		if (!this.unknownRailways.isEmpty()) {
			log.info("The following railway-types had no defaults set and were thus NOT converted:");
			for (String railwayType : this.unknownRailways) {
				log.info("- \"{}\"", railwayType);
			}
		}
		if (!this.unknownWays.isEmpty()) {
			log.info("The way-types with the following tags had no defaults set and were thus NOT converted:");
			for (String wayType : this.unknownWays) {
				log.info("- \"{}\"", wayType);
			}
		}
		log.info("= end of conversion statistics ====================");
	}

	/**
	 * Creates a MATSim link from OSM data
	 */
	protected void createLink(final Osm.Way way, final Osm.Node fromNode, final Osm.Node toNode, double length) {
		boolean oneway;
		boolean onewayReverse = false;
		double laneCapacity;
		Set<String> modes;

		// load defaults
		OsmConverterConfigGroup.OsmWayParams wayDefaultParams = getWayDefaultParams(way);
		laneCapacity = wayDefaultParams.getLaneCapacity();
		oneway = wayDefaultParams.getOneway();
		modes = new HashSet<>(wayDefaultParams.getAllowedTransportModes());

		// Overwrite defaults with OSM data
		Map<String, String> tags = way.getTags();
		String highwayValue = tags.get(Osm.Key.HIGHWAY);
		String railwayValue = tags.get(Osm.Key.RAILWAY);

		// ONEWAY
		if ("roundabout".equals(way.getTags().get(Osm.Key.JUNCTION))) {
			// if "junction" is not set in tags, get() returns null and equals() evaluates to false
			oneway = true;
		}
		String onewayTag = way.getTags().get(Osm.Key.ONEWAY);
		if (onewayTag != null) {
			switch (onewayTag) {
				case Osm.Value.YES, "true", "1" -> oneway = true;
                case "-1" -> {
					onewayReverse = true;
					oneway = false;
				}
				case "no" -> oneway = false; // may be used to overwrite defaults
			}
		}

		// FREESPEED
		double freeSpeedDefault = wayDefaultParams.getFreespeed();
		double freeSpeedForward = calculateFreeSpeed(way, true, oneway || onewayReverse, freeSpeedDefault);
		double freeSpeedBackward = calculateFreeSpeed(way, false, oneway || onewayReverse, freeSpeedDefault);

		if (config.getScaleMaxSpeed()) {
			double freeSpeedFactor = wayDefaultParams.getFreespeedFactor();
			freeSpeedForward = freeSpeedForward * freeSpeedFactor;
			freeSpeedBackward = freeSpeedBackward * freeSpeedFactor;
		}

		// LANES
		double laneCountDefault = wayDefaultParams.getLanes();

		double laneCountForward = calculateLaneCount(way, true, oneway || onewayReverse, laneCountDefault);
		Result psvLanesForward = calculateBlockingCount(way, true, oneway || onewayReverse, laneCountDefault);

		laneCountForward -= psvLanesForward.count;

		double laneCountBackward = calculateLaneCount(way, false, oneway || onewayReverse, laneCountDefault);
		Result psvLanesBackward = calculateBlockingCount(way, false, oneway || onewayReverse, laneCountDefault);

		laneCountBackward -= psvLanesBackward.count;

		// CAPACITY
		//double capacity = laneCountDefault * laneCapacity;

		// MODES
		// public transport: get relation which this way is part of, then get the relations route=* (-> the mode)
		Set<String> ptModes = new HashSet<>();
		for (Osm.Relation rel : way.getRelations().values()) {
			String osmMode = rel.getTags().get(Osm.Key.ROUTE);
			if (ptFilter.matches(rel) && osmMode != null) {
				if (osmMode.equals(Osm.Value.TROLLEYBUS)) {
					osmMode = Osm.Value.BUS;
				}
				modes.add(osmMode);
				modes.add(TransportMode.pt);
				ptModes.add(osmMode); // remember, that this is a pt mode
			}
		}

		// TURN RESTRICTIONS
		List<OsmTurnRestriction> osmTurnRestrictions = this.parseTurnRestrictions(way, modes, ptModes);


		// LENGTH
		if (length == 0.0) {
			log.warn("Attempting to create a link of length 0.0, which will mess up the routing. Fixing to 1.0!");
			length = 1.0;
		}

		if ((laneCountForward + laneCountBackward > laneCountDefault) && !oneway) {
			laneCapacity = (double) (int) laneCapacity / 2.0; // if we have more than one lane in both directions, we assume that the capacity is shared
		}

		// CREATE LINK
		// only create link, if both nodes were found, node could be null, since nodes outside a layer were dropped
		Id<Node> fromId = Id.create(fromNode.getId(), Node.class);
		Id<Node> toId = Id.create(toNode.getId(), Node.class);
		if (network.getNodes().get(fromId) != null && network.getNodes().get(toId) != null) {
			// forward link (in OSM digitization direction)
			if (!onewayReverse) {
				Id<Link> linkId = Id.create(this.id, Link.class);
				Link l = network.getFactory().createLink(linkId, network.getNodes().get(fromId), network.getNodes().get(toId));
				l.setLength(length);
				l.setFreespeed(freeSpeedForward);
				l.setCapacity(laneCountForward * laneCapacity);
				l.setNumberOfLanes(laneCountForward);
				l.setAllowedModes(modes);
				if (config.parseTurnRestrictions && !osmTurnRestrictions.isEmpty()) {
					// filter turn restrictions to those for which this link could be the from link
					List<OsmTurnRestriction> thisOsmTurnRestrictions = osmTurnRestrictions.stream()
							.filter(tr -> tr.viaNodeId() == null || tr.viaNodeId().toString().equals(toId.toString()))
							.toList();
					log.debug("Link {}: {}/{} turn restrictions attached", linkId, thisOsmTurnRestrictions.size(),
							osmTurnRestrictions.size());
					l.getAttributes().putAttribute(OSM_TURN_RESTRICTION_ATTRIBUTE_NAME, thisOsmTurnRestrictions);
				}

				network.addLink(l);
				osmIds.put(l.getId(), way.getId());
				geometryExporter.addLinkDefinition(linkId, new LinkDefinition(fromNode, toNode, way));

				// we might have dedicated lanes
				// we need to create another link for that

				if (psvLanesForward.count > 0 && modes.contains("car")) {

					Id<Link> linkIdBus = Id.create(this.id + OSM_SPECIAL_LANE, Link.class);
					Link lBus = network.getFactory().createLink(linkIdBus, network.getNodes().get(fromId), network.getNodes().get(toId));
					lBus.setLength(length);
					lBus.setFreespeed(freeSpeedForward);
					lBus.setCapacity(psvLanesForward.count * laneCapacity);
					lBus.setNumberOfLanes(psvLanesForward.count);

					// we adjust the modes allowed on other lanes
					// while this is in practice not true, it makes it better
					// for routing buses as otherwise they might end up on the 
					// other lane
					Set<String> cmodes = new HashSet<>(modes);
					if (psvLanesForward.mode.equals(Osm.Key.BUS)) {
						lBus.setAllowedModes(Set.of(Osm.Key.BUS, "pt"));
						cmodes.remove(Osm.Key.BUS);
					} else if (psvLanesForward.mode.equals(Osm.Key.TAXI)) {
						lBus.setAllowedModes(Set.of("taxi"));
						cmodes.remove(Osm.Key.TAXI);
					} else {
						lBus.setAllowedModes(Set.of(Osm.Key.BUS, "pt", Osm.Key.TAXI));
						cmodes.remove(Osm.Key.TAXI);
						cmodes.remove(Osm.Key.BUS);
					}


					l.setAllowedModes(cmodes);
					if (config.parseTurnRestrictions && !osmTurnRestrictions.isEmpty()) {
						// filter turn restrictions to those for which this link could be the from link
						List<OsmTurnRestriction> thisOsmTurnRestrictions = osmTurnRestrictions.stream()
								.filter(tr -> tr.viaNodeId() == null || tr.viaNodeId().toString().equals(toId.toString()))
								.toList();
						log.debug("Link {}: {}/{} turn restrictions attached", linkIdBus, thisOsmTurnRestrictions.size(),
								osmTurnRestrictions.size());
						lBus.getAttributes().putAttribute(OSM_TURN_RESTRICTION_ATTRIBUTE_NAME, thisOsmTurnRestrictions);
					}

					network.addLink(lBus);
					osmIds.put(lBus.getId(), way.getId());
					geometryExporter.addLinkDefinition(linkIdBus, new LinkDefinition(fromNode, toNode, way));
				}


				this.id++;
			}
			// backward link
			if (!oneway) {
				Id<Link> linkId = Id.create(this.id, Link.class);
				Link l = network.getFactory().createLink(linkId, network.getNodes().get(toId), network.getNodes().get(fromId));
				l.setLength(length);
				l.setFreespeed(freeSpeedBackward);
				l.setCapacity(laneCountBackward * laneCapacity);
				l.setNumberOfLanes(laneCountBackward);
				l.setAllowedModes(modes);
				if (config.parseTurnRestrictions) {
                    assert osmTurnRestrictions != null;
                    if (!osmTurnRestrictions.isEmpty()) {
// filter turn restrictions to those for which this link could be the from link
                        List<OsmTurnRestriction> thisOsmTurnRestrictions = osmTurnRestrictions.stream()
                                .filter(tr -> tr.viaNodeId() == null || tr.viaNodeId().toString().equals(fromId.toString()))
                                .toList();
                        log.debug("Link {}: {}/{} turn restrictions attached", linkId, thisOsmTurnRestrictions.size(),
                                osmTurnRestrictions.size());
                        l.getAttributes().putAttribute(OSM_TURN_RESTRICTION_ATTRIBUTE_NAME, thisOsmTurnRestrictions);
                    }
                }

				network.addLink(l);
				osmIds.put(l.getId(), way.getId());
				geometryExporter.addLinkDefinition(linkId, new LinkDefinition(toNode, fromNode, way));

				if (psvLanesBackward.count > 0 && modes.contains("car")) {

					Id<Link> linkIdBus = Id.create(String.valueOf(this.id) + "_spec", Link.class);
					Link lBus = network.getFactory().createLink(linkIdBus, network.getNodes().get(toId), network.getNodes().get(fromId));
					lBus.setLength(length);
					lBus.setFreespeed(freeSpeedBackward);
					lBus.setCapacity(psvLanesBackward.count * laneCapacity);
					lBus.setNumberOfLanes(psvLanesBackward.count);

					Set<String> cmodes = new HashSet<>(modes);
					if (psvLanesBackward.mode.equals(Osm.Key.BUS)) {
						lBus.setAllowedModes(Set.of(Osm.Key.BUS, "pt"));
						cmodes.remove(Osm.Key.BUS);
					} else if (psvLanesBackward.mode.equals(Osm.Key.TAXI)) {
						lBus.setAllowedModes(Set.of(Osm.Key.TAXI));
						cmodes.remove(Osm.Key.TAXI);
					} else {
						lBus.setAllowedModes(Set.of(Osm.Key.BUS, "pt", Osm.Key.TAXI));
						cmodes.remove(Osm.Key.TAXI);
						cmodes.remove(Osm.Key.BUS);
					}
					l.setAllowedModes(cmodes);
					if (config.parseTurnRestrictions && !osmTurnRestrictions.isEmpty()) {
						// filter turn restrictions to those for which this link could be the from link
						List<OsmTurnRestriction> thisOsmTurnRestrictions = osmTurnRestrictions.stream()
								.filter(tr -> tr.viaNodeId() == null || tr.viaNodeId().toString().equals(fromId.toString()))
								.toList();
						log.debug("Link {}: {}/{} turn restrictions attached", linkIdBus, thisOsmTurnRestrictions.size(),
								osmTurnRestrictions.size());
						lBus.getAttributes().putAttribute(OSM_TURN_RESTRICTION_ATTRIBUTE_NAME, thisOsmTurnRestrictions);
					}

					network.addLink(lBus);
					osmIds.put(lBus.getId(), way.getId());
					geometryExporter.addLinkDefinition(linkIdBus, new LinkDefinition(fromNode, toNode, way));
				}
				this.id++;
			}
		}
	}

	private double calculateFreeSpeed(final Osm.Way way, boolean forward, boolean isOneway, double defaultFreeSpeed) {
		double maxspeed = parseMaxspeedValueAsMs(way, Osm.Key.MAXSPEED).orElse(defaultFreeSpeed);

		// in case a specific maxspeed per direction is available this overrules the standard maxspeed
		String direction = forward ? Osm.Key.FORWARD : Osm.Key.BACKWARD;
		Optional<Double> directedMaxspeed = parseMaxspeedValueAsMs(way, Osm.Key.combinedKey(Osm.Key.MAXSPEED, direction));
		if (directedMaxspeed.isPresent()) {
			maxspeed = directedMaxspeed.get();
		}

		return maxspeed;
	}

	/**
	 * @return speed in meters per second
	 */
	private Optional<Double> parseMaxspeedValueAsMs(final Osm.Way way, String key) {
		String value = way.getTags().get(key);
		if (value == null)
			return Optional.empty();

		// take first value if more values are given
		if (value.contains(";"))
			value = value.split(";")[0];

		double conversionDivisor = 3.6;
		if (value.contains("mph")) {
			conversionDivisor = 2.237;
			value = value.replaceAll("mph", "");
		} else if (value.contains("knots")) {
			conversionDivisor = 1.944;
			value = value.replaceAll("knots", "");
		}

		if (Osm.Value.NONE.equals(value)) {
			return Optional.of(SPEED_LIMIT_NONE_KPH / conversionDivisor);
		} else if (Osm.Value.WALK.equals(value)) {
			return Optional.of(SPEED_LIMIT_WALK_KPH / conversionDivisor);
		}

		try {
			return Optional.of(Double.parseDouble(value) / conversionDivisor);
		} catch (NumberFormatException e) {
			if (!unknownMaxspeedTags.contains(value)) {
				unknownMaxspeedTags.add(value);
				log.warn("Could not parse '{}': {} (way {})", key, e.getMessage(), way.getId());
			}
			return Optional.empty();
		}
	}

	private double calculateLaneCount(final Osm.Way way, boolean forward, boolean isOneway, double defaultLaneCount) {
		double laneCount = parseLanesValue(way, Osm.Key.LANES).orElse(defaultLaneCount);

		if (!isOneway)
			laneCount /= 2;

		// in case a specific lane count per direction is available this overrules the standard lanes
		String direction = forward ? Osm.Key.FORWARD : Osm.Key.BACKWARD;
		Optional<Double> directedLaneCount = parseLanesValue(way, Osm.Key.combinedKey(Osm.Key.LANES, direction));
		if (directedLaneCount.isPresent()) {
			laneCount = directedLaneCount.get();
		}

		// sanitize
		if (laneCount < 1)
			laneCount = 1;

		return laneCount;
	}

	public record Result(double count, String mode) {
	}

	private Result calculateBlockingCount(final Osm.Way way, boolean forward, boolean isOneway, double defaultLaneCount) {
		// subtract lanes not accessible for cars
		List<String> blockingMots = Arrays.asList(Osm.Key.BUS, Osm.Key.PSV, Osm.Key.TAXI);
		double lanestoremove = 0;
		String mode = "";
		for (String blockingMot : blockingMots) {
			lanestoremove = parseReservedLanesValue(way, Osm.Key.combinedKey(Osm.Key.LANES, blockingMot)).orElse(0d);
			if (lanestoremove != 0) {
				mode = blockingMot;
				break;
			}
			lanestoremove = parseReservedLanesValue(way, Osm.Key.combinedKey(blockingMot, Osm.Key.LANES)).orElse(0d);
			if (lanestoremove != 0) {
				mode = blockingMot;
				break;
			}
		}

		// in case a specific lane count per direction is available this overrules the standard lanes
		String direction = forward ? Osm.Key.FORWARD : Osm.Key.BACKWARD;
		Optional<Double> directedLaneCount = parseLanesValue(way, Osm.Key.combinedKey(Osm.Key.LANES, direction));
		if (directedLaneCount.isPresent()) {
			lanestoremove = 0;
			for (String blockingMot : blockingMots) {
				lanestoremove = parseReservedLanesValue(way, Osm.Key.combinedKey(Osm.Key.LANES, blockingMot, direction)).orElse(0d);
				if (lanestoremove != 0) {
					mode = blockingMot;
					break;
				}
				lanestoremove = parseReservedLanesValue(way, Osm.Key.combinedKey(blockingMot, Osm.Key.LANES, direction)).orElse(0d);
				if (lanestoremove != 0) {
					mode = blockingMot;
					break;
				}
			}

		}

		if (!isOneway)
			lanestoremove /= 2;

		return new Result(lanestoremove, mode);
	}

	private Optional<Double> parseLanesValue(final Osm.Way way, String key) {
		String value = way.getTags().get(key);
		if (value == null)
			return Optional.empty();

		// it is possible that the number of blocked lanes is not reported with a number
		// but with specific designation per lane e.g., for two lanes: designated|yes
		// meaning that left lane is for pt and right one for both

		String[] parts = value.split("\\|");
		int count = 0;
		for (String p : parts) {
			// if you want to ignore leading/trailing spaces:
			if (p.trim().equals(Osm.Value.DESIGNATED) || p.trim().equals(Osm.Value.YES)) {
				count++;
			}
		}
		// if there are no lanes reserved, return 0
		if (count > 0)
			return Optional.of((double) count);

		// take first value if more values are given
		if (value.contains(";"))
			value = value.split(";")[0];

		try {
			return Optional.of(Double.parseDouble(value));
		} catch (NumberFormatException e) {
			if (!unknownLanesTags.contains(value)) {
				unknownLanesTags.add(value);
				log.warn("Could not parse '{}': {} (way {})", key, e.getMessage(), way.getId());
			}
			return Optional.empty();
		}
	}

	private Optional<Double> parseReservedLanesValue(final Osm.Way way, String key) {
		String value = way.getTags().get(key);
		if (value == null)
			return Optional.empty();

		// it is possible that the number of blocked lanes is not reported with a number
		// but with specific designation per lane e.g., for two lanes: designated|yes
		// meaning that left lane is for pt and right one for both

		String[] parts = value.split("\\|");
		int count = 0;
		for (String p : parts) {
			// if you want to ignore leading/trailing spaces:
			if (p.trim().equals(Osm.Value.DESIGNATED) || p.trim().equals(Osm.Value.YES)) {
				count++;
			}
		}
		// if there are no lanes reserved, return 0
		if (count > 0)
			return Optional.of((double) count);

		// it can happen that the value is given as a numeric
		// take first value if more values are given
		if (value.contains(";"))
			value = value.split(";")[0];

		try {
			return Optional.of(Double.parseDouble(value));
		} catch (NumberFormatException e) {
			if (!unknownLanesTags.contains(value)) {
				unknownLanesTags.add(value);
				log.warn("Could not parse '{}': {} (way {})", key, e.getMessage(), way.getId());
			}
			return Optional.empty();
		}


	}

	private void initPT() {
		ptFilter = new AllowedTagsFilter();
		ptFilter.add(Osm.ElementType.RELATION, Osm.Key.ROUTE_MASTER, Osm.Value.BUS);
		ptFilter.add(Osm.ElementType.RELATION, Osm.Key.ROUTE_MASTER, Osm.Value.TROLLEYBUS);
		ptFilter.add(Osm.ElementType.RELATION, Osm.Key.ROUTE_MASTER, Osm.Value.TRAM);
		ptFilter.add(Osm.ElementType.RELATION, Osm.Key.ROUTE_MASTER, Osm.Value.MONORAIL);
		ptFilter.add(Osm.ElementType.RELATION, Osm.Key.ROUTE_MASTER, Osm.Value.SUBWAY);
		ptFilter.add(Osm.ElementType.RELATION, Osm.Key.ROUTE_MASTER, Osm.Value.FERRY);
		ptFilter.add(Osm.ElementType.RELATION, Osm.Key.ROUTE, Osm.Value.BUS);
		ptFilter.add(Osm.ElementType.RELATION, Osm.Key.ROUTE, Osm.Value.TROLLEYBUS);
		ptFilter.add(Osm.ElementType.RELATION, Osm.Key.ROUTE, Osm.Value.RAIL);
		ptFilter.add(Osm.ElementType.RELATION, Osm.Key.ROUTE, Osm.Value.TRAIN);
		ptFilter.add(Osm.ElementType.RELATION, Osm.Key.ROUTE, Osm.Value.TRAM);
		ptFilter.add(Osm.ElementType.RELATION, Osm.Key.ROUTE, Osm.Value.LIGHT_RAIL);
		ptFilter.add(Osm.ElementType.RELATION, Osm.Key.ROUTE, Osm.Value.FUNICULAR);
		ptFilter.add(Osm.ElementType.RELATION, Osm.Key.ROUTE, Osm.Value.MONORAIL);
		ptFilter.add(Osm.ElementType.RELATION, Osm.Key.ROUTE, Osm.Value.SUBWAY);
		ptFilter.add(Osm.ElementType.WAY, Osm.Key.PSV, Osm.Value.YES);
		ptFilter.add(Osm.ElementType.WAY, Osm.Key.PSV, Osm.Value.DESIGNATED);
		ptFilter.add(Osm.ElementType.WAY, Osm.Key.BUS, Osm.Value.DESIGNATED);
		ptFilter.add(Osm.ElementType.WAY, Osm.Key.BUS, Osm.Value.DESIGNATED);

		ptDefaultParams = new OsmConverterConfigGroup.OsmWayParams("NULL", "NULL",
				1, 50 / 3.6, 1.0, 9999,
				false, Collections.singleton(TransportMode.pt));
	}

	protected boolean wayHasPublicTransit(Osm.Way way) {
		if (ptFilter.matches(way)) {
			return true;
		}
		for (Osm.Relation relation : way.getRelations().values()) {
			if (ptFilter.matches(relation)) {
				return true;
			}
		}
		return false;
	}

	protected OsmConverterConfigGroup.OsmWayParams getWayDefaultParams(Osm.Way way) {
		Map<String, String> tags = way.getTags();
		String highwayValue = tags.get(Osm.Key.HIGHWAY);
		String railwayValue = tags.get(Osm.Key.RAILWAY);

		OsmConverterConfigGroup.OsmWayParams wayDefaults = null;
		if (highwayValue != null) {
			Map<String, OsmConverterConfigGroup.OsmWayParams> highwayParams = this.wayParams.get(Osm.Key.HIGHWAY);
			if (highwayParams != null) {
				wayDefaults = highwayParams.get(highwayValue);
				if (wayDefaults == null) {
					unknownHighways.add(highwayValue);
				}
			}
		} else if (railwayValue != null) {
			Map<String, OsmConverterConfigGroup.OsmWayParams> railwayParams = this.wayParams.get(Osm.Key.RAILWAY);
			if (railwayParams != null) {
				wayDefaults = railwayParams.get(railwayValue);
				if (wayDefaults == null) {
					unknownRailways.add(railwayValue);
				}
			}
		} else {
			unknownWays.add(way.getTags().values().toString());
		}

		if (wayDefaults == null) {
			if (wayHasPublicTransit(way) && config.keepHighwaysWithPT()) {
				wayDefaults = ptDefaultParams;
			}
		}

		return wayDefaults;
	}

	/**
	 * Adds DisallowedNextLinks attributes to links. See {@link #addAttributes()}
	 * documentation as to why this cannot be done directly when creating the link.
	 */
	private void addDisallowedNextLinksAttributes() {
		network.getLinks().values().forEach(link -> {
			DisallowedNextLinks dnl = disallowedNextLinks.get(link.getId());
			if (dnl != null) {
				NetworkUtils.setDisallowedNextLinks(link, dnl);
			}
		});
	}

	/**
	 * Adds attributes to the network link. Cannot be added directly upon link creation since we need to
	 * clean the road network and attributes are not copied while filtering
	 */
	protected void addAttributes() {
		for (Link link : this.network.getLinks().values()) {
			Osm.Way way = osmData.getWays().get(osmIds.get(link.getId()));

			// way id
			link.getAttributes().putAttribute(OsmConverterConfigGroup.LINK_ATTRIBUTE_WAY_ID, Long.parseLong(way.getId().toString()));

			// all tags
			for (Map.Entry<String, String> t : way.getTags().entrySet()) {
                String key = OsmConverterConfigGroup.LINK_ATTRIBUTE_WAY_PREFIX + t.getKey();
                String val = t.getValue();
                link.getAttributes().putAttribute(key.replace("&", "AND"), val.replace("&", "AND"));
            }

			// relation info
			for (Osm.Relation rel : way.getRelations().values()) {
				// route
				String route = rel.getTags().get(Osm.Key.ROUTE);
				if (route != null) {
					String osmRouteKey = OsmConverterConfigGroup.LINK_ATTRIBUTE_RELATION_ROUTE;
					Set<String> attr = new HashSet<>(CollectionUtils.stringToSet((String) link.getAttributes().getAttribute(osmRouteKey)));
					attr.add(route);
					link.getAttributes().putAttribute(osmRouteKey, CollectionUtils.setToString(attr));
				}

				// route master
				String route_master = rel.getTags().get(Osm.Key.ROUTE_MASTER);
				if (route_master != null) {
					String osmRouteMasterKey = OsmConverterConfigGroup.LINK_ATTRIBUTE_RELATION_ROUTE_MASTER;
					Set<String> attr = new HashSet<>(CollectionUtils.stringToSet((String) link.getAttributes().getAttribute(osmRouteMasterKey)));
					attr.add(route_master);
					link.getAttributes().putAttribute(osmRouteMasterKey, CollectionUtils.setToString(attr));
				}
			}
		}
	}

	/**
	 * Makes sure that consistent routable subnetworks are created.
	 */
	protected void cleanNetwork() {
		Set<String> subnetworkModes = new HashSet<>();
		List<Network> subnetworks = new LinkedList<>();

		for (ConfigGroup params : config.getParameterSets(OsmConverterConfigGroup.RoutableSubnetworkParams.SET_NAME)) {
			OsmConverterConfigGroup.RoutableSubnetworkParams subnetworkParams = (OsmConverterConfigGroup.RoutableSubnetworkParams) params;
			final String subnetworkMode = subnetworkParams.subnetworkMode;
			final Set<String> allowedTransportModes = subnetworkParams.allowedTransportModes;
			subnetworkModes.add(subnetworkMode);
			log.info("Creating clean subnetwork for '{}' considering links of: {}", subnetworkMode, allowedTransportModes.toString());

			Network subnetwork = NetworkTools.createFilteredNetworkByLinkMode(network, allowedTransportModes);

			joinDisallowedNextLinks(subnetwork, allowedTransportModes); // if there are > 1 allowedTransportModes

			// move everything to one mode for network cleaning
			final String tmpMode = "___" + String.join("_", allowedTransportModes.toArray(new String[0])) + "___";
			copyToTmpModeAndRemoveDnlOfSubnetworkMode(subnetwork, subnetworkMode, allowedTransportModes, tmpMode);

			// clean
			DisallowedNextLinksUtils.clean(subnetwork);

			try {
				// Safe way to call TurnRestrictionsNetworkCleaner to avoid ConcurrentModificationException
				new TurnRestrictionsNetworkCleaner().run(subnetwork, tmpMode);
			} catch (Exception e) {
				log.warn("Error during turn restrictions network cleaning for mode {}: {}", tmpMode, e.getMessage());
				log.warn("Continuing without turn restrictions cleaning for this mode");
			}

			// remove all links without tmpMode
			subnetwork.getLinks().values().stream()
					.filter(link -> !link.getAllowedModes().contains(tmpMode))
					.map(Link::getId)
					.toList()
					.forEach(subnetwork::removeLink);
			NetworkUtils.removeLinksWithoutModes(subnetwork);
			NetworkUtils.removeNodesWithoutLinks(subnetwork);

			// assign subnetworkMode to remaining links of tmpMode
			DisallowedNextLinksUtils.copy(subnetwork, tmpMode, subnetworkMode);
			Set<String> subnetworkModeSingleton = Collections.singleton(subnetworkMode);
			subnetwork.getLinks().values().stream()
					.filter(link -> link.getAllowedModes().contains(tmpMode)) // should not exclude anything anymore
					.forEach(link -> link.setAllowedModes(subnetworkModeSingleton));
			DisallowedNextLinksUtils.clean(subnetwork);

			subnetworks.add(subnetwork);
		}

		Set<String> remainingModes = new HashSet<>();
		for (Link link : network.getLinks().values()) {
			remainingModes.addAll(link.getAllowedModes());
		}
		remainingModes.removeAll(subnetworkModes);

		log.info("Creating remaining network with modes: {}", remainingModes.toString());
		Network remainingNetwork = NetworkTools.createFilteredNetworkByLinkMode(network, remainingModes);

		for (Link link : remainingNetwork.getLinks().values()) {
			Set<String> newAllowedModes = new HashSet<>(remainingModes);
			newAllowedModes.retainAll(link.getAllowedModes());
			link.setAllowedModes(newAllowedModes);
		}

		subnetworks.add(remainingNetwork);

		log.info("Creating combined network");
		Network combinedNetwork = NetworkUtils.createNetwork();
		subnetworks.forEach(n -> NetworkTools.integrateNetwork(combinedNetwork, n, true));
		DisallowedNextLinksUtils.clean(combinedNetwork);

		this.network = combinedNetwork;
	}

	/**
	 * Safe method to run turn restrictions network cleaner to avoid ConcurrentModificationException
	 */
	private void safeTurnRestrictionsNetworkCleaner(Network network, String mode) {
		// First make a defensive copy of all the links that have the mode
		Set<Link> linksWithMode = network.getLinks().values().stream()
			.filter(link -> link.getAllowedModes().contains(mode))
			.collect(Collectors.toSet());

		// Process each link safely
		for (Link link : linksWithMode) {
			if (!network.getLinks().containsKey(link.getId())) {
				continue; // Skip if link was removed in previous iterations
			}

			// Create a cleaner for each link to avoid concurrent modification
			try {
				DisallowedNextLinks dnl = NetworkUtils.getDisallowedNextLinks(link);
				if (dnl != null && !dnl.getDisallowedLinkSequences(mode).isEmpty()) {
					// Handle each link's disallowed sequences individually
					DisallowedNextLinksUtils.clean(network);
				}
			} catch (Exception e) {
				log.warn("Error processing turn restrictions for link {}: {}", link.getId(), e.getMessage());
			}
		}
	}

	// Turn Restrictions

	@Nullable
	private List<OsmTurnRestriction> parseTurnRestrictions(final Osm.Way way, Set<String> modes, Set<String> ptModes) {

		if (!config.parseTurnRestrictions) {
			return null;
		}

		List<OsmTurnRestriction> osmTurnRestrictions = new ArrayList<>();
		for (Osm.Relation relation : way.getRelations().values()) {

			Map<String, String> relationTags = relation.getTags();

			// we only consider this relation, if
			// - it is a turn restriction relation and
			// - this way is the "from" link
			if (!(Osm.Key.RESTRICTION.equals(relationTags.get(Osm.Key.TYPE))
					&& relation.getMemberRoles(way).contains(Osm.Value.FROM))) {
				continue;
			}

			// identify modes
			Set<String> restrictionModes = new HashSet<>(modes);
			// remove except modes
			String exceptOsmModesString = relationTags.get(Osm.Key.EXCEPT);
			if (exceptOsmModesString != null) {
				Set<String> exceptOsmModes = Set.of(exceptOsmModesString.split(";"));
				Set<String> exceptModes = exceptOsmModes.stream()
						.map(m -> OSM_2_MATSIM_MODE_MAP.getOrDefault(m, m))
						.collect(Collectors.toSet());
				restrictionModes.removeAll(exceptModes);
				// remove pt, if all pt modes are excluded
				if (!ptModes.isEmpty() && exceptModes.containsAll(ptModes)) {
					restrictionModes.remove(TransportMode.pt);
					log.info("OSM:{} Removed pt from restricted modes: {}", way.getId(), ptModes);
				}
			}

			// identify restriction type and eventually add modes
			OsmTurnRestriction.RestrictionType restrictionType = null;
			for (String suffix : TURN_RESTRICTION_KEY_SUFFIXES) {
				String restrictionTypeString = relationTags.get(Osm.Key.RESTRICTION + suffix);
				if (restrictionTypeString != null) {

					// add restriction type
					if (restrictionTypeString.startsWith(Osm.Key.PROHIBITORY_RESTRICTION_PREFIX)) {
						restrictionType = OsmTurnRestriction.RestrictionType.PROHIBITIVE;
					} else if (restrictionTypeString.startsWith(Osm.Key.MANDATORY_RESTRICTION_PREFIX)) {
						restrictionType = OsmTurnRestriction.RestrictionType.MANDATORY;
					}

					// add explicit modes, if
					// - suffix specified it and
					// - it is a MATSim mode
					if (suffix.length() > 1) {
						String osmMode = suffix.substring(1);
						String mode = OSM_2_MATSIM_MODE_MAP.get(osmMode);
						if (mode == null) {
							// skip this, if not one of MATSim modes
							restrictionType = null;
							continue;
						}
						restrictionModes.add(mode);
					}

					break; // take first one
				}
			}
			if (restrictionType == null) {
				log.warn("Could not identify turn restriction relation: https://www.openstreetmap.org/relation/{}",
						relation.getId());
				continue;
			}

			// create intermediate turn restriction record
			List<Id<Osm.Way>> nextWayIds = new ArrayList<>();
			Id<Osm.Way> toWayId = null;
			Id<Osm.Node> viaNodeId = null;
			for (Osm.Element element : relation.getMembers()) {
				List<String> memberRoles = relation.getMemberRoles(element);
				if (element instanceof Osm.Way wayElement) {
					if (memberRoles.contains(Osm.Value.TO)) {
						toWayId = wayElement.getId();
					} else if (memberRoles.contains(Osm.Value.VIA)) {
						nextWayIds.add(wayElement.getId());
					}
				} else if (element instanceof Osm.Node nodeElement && memberRoles.contains(Osm.Value.VIA)) {
					viaNodeId = nodeElement.getId();
				}
			}
			nextWayIds.add(toWayId);
			osmTurnRestrictions.add(new OsmTurnRestriction(restrictionModes, nextWayIds, viaNodeId, restrictionType));
		}

		if (osmTurnRestrictions.isEmpty()) {
			return Collections.emptyList();
		}

		return osmTurnRestrictions;
	}

	private void attachTurnRestrictionsAsDisallowedNextLinks() {

		if (!config.parseTurnRestrictions) {
			return;
		}

		Map<Id<Link>, DisallowedNextLinks> map = network.getLinks().entrySet().parallelStream().map(e -> {
					Link link = e.getValue();

					// Note, that these turn restriction attributes are attached to both forward and
					// reverse MATSim links, but only one is valid.
					@SuppressWarnings("unchecked")
					List<OsmTurnRestriction> osmTurnRestrictions = (List<OsmTurnRestriction>) link.getAttributes()
							.getAttribute(OSM_TURN_RESTRICTION_ATTRIBUTE_NAME);
					if (osmTurnRestrictions == null) {
						return null;
					}

					// create DisallowedNextLinks
					DisallowedNextLinks dnl = null;
					for (OsmTurnRestriction tr : osmTurnRestrictions) {

						// find next link ids from next way ids
						List<Id<Link>> nextLinkIds = findLinkIds(this.wayLinkMap, this.network, link.getToNode(),
								tr.nextWayIds);
						if (nextLinkIds.isEmpty()) {
							continue;
						}

						// find link id lists to disallow
						List<List<Id<Link>>> disallowedNextLinkIdLists = new ArrayList<>();
						if (tr.restrictionType.equals(OsmTurnRestriction.RestrictionType.PROHIBITIVE)) {
							disallowedNextLinkIdLists.add(nextLinkIds);
						} else if (tr.restrictionType.equals(OsmTurnRestriction.RestrictionType.MANDATORY)) {
							// we need to exclude all other links originating from fromWay's toNode
							link.getToNode().getOutLinks().values().stream()
									.map(Link::getId)
									.filter(lId -> !lId.equals(nextLinkIds.getFirst()))
									.forEach(lId -> disallowedNextLinkIdLists.add(List.of(lId)));
						}

						// attach DisallowedNextLinks objects
						if (dnl == null && !disallowedNextLinkIdLists.isEmpty()) {
							dnl = new DisallowedNextLinks();
							log.debug("Link {}: modes={} disallowedNextLinkIdLists={}", e.getKey(),
									Arrays.toString(tr.modes.toArray()), Arrays.toString(disallowedNextLinkIdLists.toArray()));
						}
						for (List<Id<Link>> disallowedNextLinkIds : disallowedNextLinkIdLists) {
							for (String mode : tr.modes) {
								dnl.addDisallowedLinkSequence(mode, disallowedNextLinkIds);
							}
						}
					}

					// remove attribute
					link.getAttributes().removeAttribute(OSM_TURN_RESTRICTION_ATTRIBUTE_NAME);

					if (dnl != null) {
						return Map.entry(e.getKey(), dnl);
					}
					return null;
				})
				.filter(Objects::nonNull)
				.collect(Collectors.toMap(Entry::getKey, Entry::getValue));

		this.disallowedNextLinks.putAll(map);
	}

	// Statics

	private static void copyToTmpModeAndRemoveDnlOfSubnetworkMode(Network subnetwork, String subnetworkMode,
																  Set<String> allowedTransportModes, String tmpMode) {
		Verify.verify(subnetwork.getLinks().values().stream()
				.noneMatch(link -> link.getAllowedModes().contains(tmpMode)));
		subnetwork.getLinks().values().forEach(link -> NetworkUtils.addAllowedMode(link, tmpMode));
		// copy remaining DNLs from allowed modes to tmpMode
		for (String allowedMode : allowedTransportModes) {
			DisallowedNextLinksUtils.copy(subnetwork, allowedMode, tmpMode);
			// remove DNLs from subnetworkMode, as final mode will be subnetworkMode
			for (Link link : subnetwork.getLinks().values()) {
				DisallowedNextLinks dnl = NetworkUtils.getDisallowedNextLinks(link);
				if (dnl != null) {
					dnl.removeDisallowedLinkSequences(subnetworkMode);
					if (dnl.isEmpty()) {
						NetworkUtils.removeDisallowedNextLinks(link);
					}
				}
			}
		}
	}

	/**
	 * If we want to create a routable subnetwork on multiple modes, adding all turn
	 * restrictions of both modes might be too restrictive on the joined network, so
	 * we remove turn restrictions that would have previously forbidden to transfer
	 * between links that have different allowed modes (only considering
	 * allowedTransportModes).
	 *
	 * @param subnetwork
	 * @param allowedTransportModes
	 */
	private static void joinDisallowedNextLinks(Network subnetwork, Set<String> allowedTransportModes) {

		if (allowedTransportModes.size() > 1) {
			return; // skip, as there are no DNLs to join from multiple modes
		}

		// Note: DNLs are not cleaned yet!
		// When joining subnetworks from allowedTransportModes, remove all DNL link
		// sequences, that would prohibit traversing between links with *different*
		// allowedModes.
		List<List<Id<Link>>> linkConnectionsToRemove = new ArrayList<>();
		for (Link link : subnetwork.getLinks().values()) {
			Set<String> modes = new HashSet<>(link.getAllowedModes());
			modes.retainAll(allowedTransportModes);

			for (Link otherLink : link.getToNode().getOutLinks().values()) {
				Set<String> otherModes = new HashSet<>(otherLink.getAllowedModes());
				otherModes.retainAll(allowedTransportModes);

				if (!modes.equals(otherModes)) {
					linkConnectionsToRemove.add(List.of(link.getId(), otherLink.getId()));
				}
			}
		}
		// remove some link connections from DNLs
		// a) link connection starts from link 1 -> link2 is first in a DNL sequence
		linkConnectionsToRemove.forEach(t -> {
			Link link1 = subnetwork.getLinks().get(t.get(0));
			Link link2 = subnetwork.getLinks().get(t.get(1));

			DisallowedNextLinks dnl = NetworkUtils.getDisallowedNextLinks(link1);
			if (dnl != null) {
				for (String allowedMode : allowedTransportModes) {
					List<List<Id<Link>>> linkSequences = new ArrayList<>(
							dnl.getDisallowedLinkSequences(allowedMode));
					boolean removedLinkSequence = false;
					for (ListIterator<List<Id<Link>>> it = linkSequences.listIterator(); it.hasNext(); ) {
						List<Id<Link>> linkSequence = it.next();
						if (linkSequence.size() == 1 && linkSequence.getFirst().equals(link2.getId())) {
							it.remove();
							removedLinkSequence = true;
							log.warn("Removed link sequence {} from {} for {}", linkSequence.toString(),
									link1.getId(), allowedMode); // ! DEBUG
						}
					}
					if (removedLinkSequence) {
						dnl.removeDisallowedLinkSequences(allowedMode);
						linkSequences.forEach(ls -> dnl.addDisallowedLinkSequence(allowedMode, ls));
					}
				}
			}
		});
		// b) both links of link connection are in the link sequence of another DNL
		subnetwork.getLinks().values().parallelStream().forEach(link -> {
			DisallowedNextLinks dnl = NetworkUtils.getDisallowedNextLinks(link);
			if (dnl != null) {
				for (String allowedMode : allowedTransportModes) {
					List<List<Id<Link>>> linkSequences = new ArrayList<>(
							dnl.getDisallowedLinkSequences(allowedMode));
					boolean removedLinkSequence = false;
					for (ListIterator<List<Id<Link>>> it = linkSequences.listIterator(); it.hasNext(); ) {
						List<Id<Link>> linkSequence = it.next();
						for (List<Id<Link>> linkConnection : linkConnectionsToRemove) { // not very efficient
							if (Collections.indexOfSubList(linkSequence, linkConnection) >= 0) {
								it.remove();
								removedLinkSequence = true;
								log.warn("Removed link sequence {} from {} for {}", linkSequence.toString(),
										link.getId(), allowedMode); // ! DEBUG
							}
						}
					}
					if (removedLinkSequence) {
						dnl.removeDisallowedLinkSequences(allowedMode);
						linkSequences.forEach(
								linkSequence -> dnl.addDisallowedLinkSequence(allowedMode, linkSequence));
					}
				}
			}
		});
	}

	/**
	 * Finds list of link ids starting from {@code lastNode} from list of OSM way
	 * ids.
	 *
	 * @param wayLinkMap
	 * @param network
	 * @param lastNode
	 * @param wayIds
	 * @return
	 */
	protected static List<Id<Link>> findLinkIds(Map<Id<Osm.Way>, List<Id<Link>>> wayLinkMap, Network network,
												Node lastNode, List<Id<Osm.Way>> wayIds) {

		Map<Id<Link>, ? extends Link> links = network.getLinks();
		List<Entry<Id<Osm.Way>, Id<Link>>> linkIds = new ArrayList<>();

		log.debug("Ways {}: Checking for valid turn restriction", Arrays.toString(wayIds.toArray()));

		for (Id<Osm.Way> wayId : wayIds) {

			// for every link id, that could stem from this way
			List<Id<Link>> linkIdCandidates = wayLinkMap.get(wayId);
			if (linkIdCandidates == null) {
				// requested way id has no link ids -> turn restriction is incomplete
				log.debug("Invalid: Way {} has no links in MATSim network.", wayId);
				return Collections.emptyList();
			} else {
				linkIdCandidates = new ArrayList<>(linkIdCandidates);
			}

			// Ensure, there is only one link candidate starting from lastNode. If there are
			// two, the turn restriction is not valid as we do not know which one to take as
			// the start.
			int linksStartingAtLastNode = 0;
			for (Id<Link> linkIdCandidate : linkIdCandidates) {
				Link linkCandidate = links.get(linkIdCandidate);
				if (lastNode.getId().equals(linkCandidate.getFromNode().getId())) {
					linksStartingAtLastNode++;
				}
				if (linksStartingAtLastNode > 1) {
					log.debug("Invalid: Way {} has multiple links starting from lastNode.", wayId);
					return Collections.emptyList();
				}
			}


			log.debug("  Way {} might belong to valid turn restriction with links {}", wayId, Arrays.toString(linkIdCandidates.toArray()));

			// find link sequence from way
			Link nextLink;
			do {
				nextLink = null;

				// find next link id
				for (Iterator<Id<Link>> linkIdIt = linkIdCandidates.listIterator(); linkIdIt.hasNext(); ) {
					Id<Link> linkIdCandidate = linkIdIt.next();
					Link linkCandidate = links.get(linkIdCandidate);

					if (lastNode.getId().equals(linkCandidate.getFromNode().getId())) {
						nextLink = linkCandidate; // subsequent link found from lastNode
						linkIds.add(Map.entry(wayId, linkIdCandidate));
						log.debug("  Way {}: Next link {}", wayId, linkIdCandidate);

						// remove this link from candidates
						linkIdIt.remove();

						lastNode = linkCandidate.getToNode();
						break; // we need to remove the reverse link first, before continuing
					}
				}

				// remove reverse link of found link, so it is not found in next iteration
				if (nextLink != null) {
					for (Iterator<Id<Link>> linkIdIt = linkIdCandidates.listIterator(); linkIdIt.hasNext(); ) {
						Id<Link> linkIdCandidate = linkIdIt.next();
						Link linkCandidate = links.get(linkIdCandidate);

						if (isReverse(nextLink, linkCandidate)) {
							linkIdIt.remove();
							log.debug("  Way {}: Next link {} -> removed reverse link {}", wayId, nextLink.getId(),
									linkIdCandidate);
						}
					}
				}

				// repeat until no new next link is found
				log.debug("  Way {}: Next link {} -> remaining linkIdCandidates={}", wayId,
						nextLink == null ? "(null)" : nextLink.getId(),
						Arrays.toString(linkIdCandidates.toArray()));
			} while (nextLink != null && !linkIdCandidates.isEmpty());

		}

		// turn restrictions are only valid, if at least one link is found for each way
		if (linkIds.size() < wayIds.size()) {
			log.debug("Invalid: Ways {} did not have at least one MATSim link each.",
					Arrays.toString(wayIds.toArray()));
			return Collections.emptyList(); // turn restriction could not be matched to MATSim network
		}

		log.debug("Valid: Ways {} have valid turn restriction.", Arrays.toString(wayIds.toArray()));
		linkIds.forEach(e -> log.debug("  Way {}: {}", e.getKey(), e.getValue()));

		return linkIds.stream().map(Entry::getValue).toList();
	}

	private static boolean isReverse(Link f, Link r) {
		return f.getToNode().getId().equals(r.getFromNode().getId())
				&& f.getFromNode().getId().equals(r.getToNode().getId());
	}

}

