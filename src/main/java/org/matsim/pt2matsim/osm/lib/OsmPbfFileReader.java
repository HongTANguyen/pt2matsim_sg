/* *********************************************************************** *
 * project: org.matsim.*
 *                                                                         *
 * *********************************************************************** *
 *                                                                         *
 * copyright       : (C) 2013 by the members listed in the COPYING,        *
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
 * *********************************************************************** */

package org.matsim.pt2matsim.osm.lib;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.matsim.api.core.v01.Coord;
import org.matsim.core.utils.misc.Counter;

import com.slimjars.dist.gnu.trove.iterator.TLongObjectIterator;
import de.topobyte.osm4j.core.dataset.InMemoryMapDataSet;
import de.topobyte.osm4j.core.dataset.MapDataSetLoader;
import de.topobyte.osm4j.core.model.iface.*;
import de.topobyte.osm4j.pbf.seq.PbfIterator;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.UncheckedIOException;
import java.util.HashMap;
import java.util.Map;

import static org.matsim.pt2matsim.osm.lib.Osm.ElementType.*;

/**
 * Reads OpenStreetMap data from a PBF file and converts it to the MATSim {@link OsmData} format.
 *
 * @author polettif (based on OsmXmlFileReader)
 */
public class OsmPbfFileReader {

    private static final Logger log = LogManager.getLogger(OsmPbfFileReader.class);
    private final OsmData osmData;

    /**
     * Constructor.
     *
     * @param osmData the container to store parsed OSM data in
     */
    public OsmPbfFileReader(OsmData osmData) {
        this.osmData = osmData;
    }

    /**
     * Parse the PBF file and load OSM data into the OsmData container.
     *
     * @param filename path to the PBF file
     */
    public void readFile(String filename) {
        log.info("Loading OSM data from PBF file {}...", filename);

        final Counter nodeCounter = new Counter("node ");
        final Counter wayCounter = new Counter("way ");
        final Counter relationCounter = new Counter("relation ");

        try (InputStream inputStream = new FileInputStream(filename)) {
            // Create PBF iterator and load data into memory
            PbfIterator pbfIterator = new PbfIterator(inputStream, true);
            InMemoryMapDataSet mapDataSet = MapDataSetLoader.read(pbfIterator, true, true, true);

            log.info("Processing nodes...");
            // Process nodes
            TLongObjectIterator<OsmNode> nodeIt = mapDataSet.getNodes().iterator();
            while (nodeIt.hasNext()) {
                nodeIt.advance();
                OsmNode osmNode = nodeIt.value();
                nodeCounter.incCounter();

                long id = osmNode.getId();
                // Create parsed node object
                OsmXmlFileReader.ParsedNode node = new OsmXmlFileReader.ParsedNode(
                    id,
                    new Coord(osmNode.getLongitude(), osmNode.getLatitude())
                );

                // Add tags
                if (osmNode.getNumberOfTags() > 0) {
                    for (int i = 0; i < osmNode.getNumberOfTags(); i++) {
                        node.tags.put(
                            StringCache.get(osmNode.getTag(i).getKey()),
                            StringCache.get(osmNode.getTag(i).getValue())
                        );
                    }
                }

                // Pass to OSM data container
                osmData.handleParsedNode(node);
            }
            nodeCounter.printCounter();

            log.info("Processing ways...");
            // Process ways
            TLongObjectIterator<OsmWay> wayIt = mapDataSet.getWays().iterator();
            while (wayIt.hasNext()) {
                wayIt.advance();
                OsmWay osmWay = wayIt.value();
                wayCounter.incCounter();

                long id = osmWay.getId();
                // Create parsed way object
                OsmXmlFileReader.ParsedWay way = new OsmXmlFileReader.ParsedWay(id);

                // Add nodes to way
                for (int i = 0; i < osmWay.getNumberOfNodes(); i++) {
                    way.nodes.add(osmWay.getNodeId(i));
                }

                // Add tags
                if (osmWay.getNumberOfTags() > 0) {
                    for (int i = 0; i < osmWay.getNumberOfTags(); i++) {
                        way.tags.put(
                            StringCache.get(osmWay.getTag(i).getKey()),
                            StringCache.get(osmWay.getTag(i).getValue())
                        );
                    }
                }

                // Pass to OSM data container
                osmData.handleParsedWay(way);
            }
            wayCounter.printCounter();

            log.info("Processing relations...");
            // Process relations
            TLongObjectIterator<OsmRelation> relationIt = mapDataSet.getRelations().iterator();
            while (relationIt.hasNext()) {
                relationIt.advance();
                OsmRelation osmRelation = relationIt.value();
                relationCounter.incCounter();

                long id = osmRelation.getId();
                // Create parsed relation object
                OsmXmlFileReader.ParsedRelation relation = new OsmXmlFileReader.ParsedRelation(id);

                // Add members
                for (int i = 0; i < osmRelation.getNumberOfMembers(); i++) {
                    OsmRelationMember member = osmRelation.getMember(i);

                    Osm.ElementType type = switch (member.getType()) {
                        case Node -> NODE;
                        case Way -> WAY;
                        case Relation -> RELATION;
                    };

                    relation.members.add(new OsmXmlFileReader.ParsedRelationMember(
                        type,
                        member.getId(),
                        StringCache.get(member.getRole())
                    ));
                }

                // Add tags
                if (osmRelation.getNumberOfTags() > 0) {
                    for (int i = 0; i < osmRelation.getNumberOfTags(); i++) {
                        relation.tags.put(
                            StringCache.get(osmRelation.getTag(i).getKey()),
                            StringCache.get(osmRelation.getTag(i).getValue())
                        );
                    }
                }

                // Pass to OSM data container
                osmData.handleParsedRelation(relation);
            }
            relationCounter.printCounter();

            // Finalize OSM data
            log.info("Building OSM network graph...");
            osmData.buildMap();

            log.info("Finished loading OSM data: {} nodes, {} ways, {} relations", nodeCounter.getCounter(), wayCounter.getCounter(), relationCounter.getCounter());

        } catch (IOException e) {
            throw new UncheckedIOException("Error reading OSM PBF file: " + filename, e);
        }
    }

    /**
     * String cache utility for efficient memory usage.
     */
    public static class StringCache {
        private static final Map<String, String> CACHE = new HashMap<>();

        public static String get(String s) {
            if (s == null) return null;
            String cached = CACHE.get(s);
            if (cached == null) {
                CACHE.put(s, s);
                return s;
            }
            return cached;
        }
    }
}
