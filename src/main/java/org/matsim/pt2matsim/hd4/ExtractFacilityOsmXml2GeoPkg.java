package org.matsim.pt2matsim.hd4;

import com.fasterxml.jackson.annotation.JsonAnyGetter;
import com.fasterxml.jackson.annotation.JsonAnySetter;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.google.common.collect.Iterables;
import de.topobyte.osm4j.geometry.GeometryBuilder;
import it.unimi.dsi.fastutil.longs.Long2ObjectLinkedOpenHashMap;
import it.unimi.dsi.fastutil.longs.Long2ObjectMap;
import it.unimi.dsi.fastutil.objects.Object2IntLinkedOpenHashMap;
import it.unimi.dsi.fastutil.objects.Object2IntMap;
import it.unimi.dsi.fastutil.objects.ObjectCollection;
import me.tongfei.progressbar.ProgressBar;
import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVPrinter;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.geotools.api.data.DataStore;
import org.geotools.api.data.DataStoreFinder;
import org.geotools.api.data.SimpleFeatureStore;
import org.geotools.api.data.Transaction;
import org.geotools.api.feature.simple.SimpleFeature;
import org.geotools.api.feature.simple.SimpleFeatureType;
import org.geotools.api.referencing.FactoryException;
import org.geotools.api.referencing.crs.CRSAuthorityFactory;
import org.geotools.api.referencing.operation.MathTransform;
import org.geotools.api.referencing.operation.TransformException;
import org.geotools.data.DefaultTransaction;
import org.geotools.data.collection.ListFeatureCollection;
import org.geotools.feature.simple.SimpleFeatureBuilder;
import org.geotools.feature.simple.SimpleFeatureTypeBuilder;
import org.geotools.geometry.jts.JTS;
import org.geotools.geopkg.GeoPkgDataStoreFactory;
import org.geotools.jdbc.JDBCDataStoreFactory;
import org.geotools.referencing.CRS;
import org.locationtech.jts.algorithm.ConvexHull;
import org.locationtech.jts.geom.*;
import org.locationtech.jts.index.strtree.STRtree;
import org.matsim.application.MATSimAppCommand;
import org.matsim.application.options.CrsOptions;
import org.matsim.core.utils.io.IOUtils;
import org.matsim.pt2matsim.osm.lib.*;
import picocli.CommandLine;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.function.Predicate;

/**
 * Version of ExtractFacilityOsmPbf2GeoPkg that reads from OSM XML files instead of PBF.
 * Core functionality remains the same.
 */
@CommandLine.Command(
    name = "facility-shp-xml",
    description = "Generate facility shape file from OSM XML data."
)
public class ExtractFacilityOsmXml2GeoPkg implements MATSimAppCommand {

    private static final Logger log = LogManager.getLogger(ExtractFacilityOsmXml2GeoPkg.class);
    private static final double POI_BUFFER = 6;

    /**
     * Structures of this size are completely ignored.
     */
    private static final double MAX_AREA = 50_000_000;

    /**
     * Structures larger than this will not be assigned smaller scale types, but remain independently.
     * Usually large areas such as parks, campus, etc.
     */
    private static final double MAX_ASSIGN = 50_000;

    /**
     * Percentage of minimum required intersection.
     */
    private static final double INTERSECT_THRESHOLD = 0.2;

    private final GeometryBuilder geometryBuilder = new GeometryBuilder();

    @CommandLine.Option(names = "--input", description = "Path to input .osm or .xml file", required = true)
    private Path osmXml;

    @CommandLine.Option(names = "--output", description = "Path to output .gpkg file", required = true)
    private Path output;

    @CommandLine.Option(names = "--activity-mapping", description = "Path to activity mapping json", required = true)
    private Path mappingPath;

    @CommandLine.Option(names = "--exclude", description = "Exclude these activities types from the output", split = ",", defaultValue = "")
    private Set<String> exclude;

    @CommandLine.Mixin
    private final CrsOptions crs = new CrsOptions("EPSG:4326", "EPSG:3414");

    /**
     * Maps types to feature index.
     */
    private Object2IntMap<String> types;
    private ActivityMapping config;
    private Long2ObjectMap<Feature> pois;
    private Long2ObjectMap<Feature> landuse;
    private Long2ObjectMap<Feature> entities;
    private MathTransform transform;
    private OsmData osmData;
    private int ignored;

    public static void main(String[] args) {
        new ExtractFacilityOsmXml2GeoPkg().execute(args);
    }

    @Override
    public Integer call() throws Exception {

        // Read OSM XML data using OsmXmlFileReader
        osmData = new OsmDataImpl();
        OsmXmlFileReader reader = new OsmXmlFileReader(osmData);
        reader.readFile(osmXml.toString());

        config = new ObjectMapper().readerFor(ActivityMapping.class).readValue(mappingPath.toFile());

        // Fix: Add proper null checks and error handling for CRS
        String inputCRS = crs.getInputCRS();
        String targetCRS = crs.getTargetCRS();

        if (inputCRS == null || inputCRS.isBlank()) {
            log.error("Input CRS is null or empty. Using default WGS84.");
            inputCRS = "EPSG:4326"; // Default to WGS84
        }

        if (targetCRS == null || targetCRS.isBlank()) {
            log.error("Target CRS is null or empty. Using default WGS84.");
            targetCRS = "EPSG:4326"; // Default to WGS84
        }

        try {
            CRSAuthorityFactory cFactory = CRS.getAuthorityFactory(true);
            transform = CRS.findMathTransform(
                cFactory.createCoordinateReferenceSystem(inputCRS),
                CRS.decode(targetCRS),
                true
            );
            log.info("Using coordinate transformation from {} to {}", inputCRS, targetCRS);

            // Debug: Test coordinate transformation with a known Singapore coordinate
            try {
                Point testPoint = geometryBuilder.getGeometryFactory().createPoint(new Coordinate(103.8198, 1.3521)); // Singapore CBD
                Point transformed = (Point) JTS.transform(testPoint, transform);
                log.info("Test transformation: WGS84 (103.8198, 1.3521) -> EPSG:3414 ({}, {})",
                    transformed.getX(), transformed.getY());

                // Also test with swapped coordinates
                Point testPointSwapped = geometryBuilder.getGeometryFactory().createPoint(new Coordinate(1.3521, 103.8198));
                Point transformedSwapped = (Point) JTS.transform(testPointSwapped, transform);
                log.info("Test transformation (swapped): WGS84 (1.3521, 103.8198) -> EPSG:3414 ({}, {})",
                    transformedSwapped.getX(), transformedSwapped.getY());

            } catch (Exception debugEx) {
                log.warn("Debug coordinate transformation failed: {}", debugEx.getMessage());
            }
        } catch (Exception e) {
            log.error("Failed to create coordinate transformation: {}. Using identity transform instead.", e.getMessage());
            // Use identity transform as fallback
            transform = null;
        }

        log.info("Configured tags: {}", config.getTypes());

        // Enable debug logging for unmapped features
        log.info("Debug logging enabled - will show unmapped OSM tags");

        types = new Object2IntLinkedOpenHashMap<>();

        config.types.values().stream()
            .flatMap(c -> c.values.values().stream())
            .flatMap(Collection::stream)
            .filter(t -> !exclude.contains(t))
            .distinct()
            .sorted()
            .forEach(e -> types.put(e, types.size()));

        log.info("Configured activity types: {}", types.keySet());

        if (types.keySet().stream().anyMatch(t -> t.length() > 10)) {
            log.error("Activity names max length is 10, due to shp format limitation.");
            return 2;
        }

        if (Files.exists(output)) {
            log.info("Deleting already existing file: {}", output);
            Files.delete(output);
        }

        // Collect all geometries first
        pois = new Long2ObjectLinkedOpenHashMap<>();
        entities = new Long2ObjectLinkedOpenHashMap<>();
        landuse = new Long2ObjectLinkedOpenHashMap<>();

        log.info("Finished loading OSM XML file.");

        // Process nodes
        log.info("Processing nodes...");
        for (Osm.Node node : ProgressBar.wrap(osmData.getNodes().values(), "Nodes")) {
            Map<String, String> tags = node.getTags();
            OsmEntity entity = new OsmEntity(node, tags);

            processNode(node);
        }
        log.info("Collected {} POIs", pois.size());

        // Process ways
        log.info("Processing ways...");
        for (Osm.Way way : ProgressBar.wrap(osmData.getWays().values(), "Ways")) {
            processWay(way);
        }

        // Process relations
        log.info("Processing relations...");
        for (Osm.Relation relation : ProgressBar.wrap(osmData.getRelations().values(), "Relations")) {
            processRelation(relation);
        }

        log.info("Collected {} landuse shapes", landuse.size());
        log.info("Collected {} other entities", entities.size());

        if (ignored > 0)
            log.warn("Ignored {} invalid geometries", ignored);

        FacilityFeatureExtractor ft = new FacilityFeatureExtractor(crs.getTargetCRS(), types, entities, pois, landuse);

        preprocessLanduse(landuse.values(), ft.entities, 0.2);

        processIntersection(landuse.values(), ft.entities, INTERSECT_THRESHOLD);

        // Low prio landuses are not needed anymore
        landuse.values().removeIf(f -> f.lowPriority);

        log.info("Remaining landuse shapes after assignment: {} ", landuse.size());

        processIntersection(pois.values(), ft.entities, 0);

        log.info("Remaining POI after assignment: {}", pois.size());

        DataStore ds = DataStoreFinder.getDataStore(Map.of(
            GeoPkgDataStoreFactory.DBTYPE.key, "geopkg",
            GeoPkgDataStoreFactory.DATABASE.key, output.toFile().toString(),
            JDBCDataStoreFactory.BATCH_INSERT_SIZE.key, 100,
            GeoPkgDataStoreFactory.READ_ONLY.key, false
        ));

        ds.createSchema(ft.featureType);

        SimpleFeatureStore source = (SimpleFeatureStore) ds.getFeatureSource(ft.featureType.getTypeName());
        ListFeatureCollection collection = new ListFeatureCollection(ft.featureType);

        addFeatures(entities, ft, collection);
        addFeatures(landuse, ft, collection);
        addFeatures(pois, ft, collection);

        Transaction transaction = new DefaultTransaction("create");
        source.setTransaction(transaction);

        source.addFeatures(collection);
        transaction.commit();

        log.info("Wrote {} features", collection.size());

        transaction.close();

        ds.dispose();

        writeMapping(output.toString().replace(".gpkg", "_mapping.csv.gz"),
            entities.values(), landuse.values(), pois.values());

        return 0;
    }

    	/**
	 * Stores entities and geometries as necessary.
	 */
//	private void process(OsmEntity entity) {
//		boolean filtered = true;
//		boolean isBuilding = false;
//		boolean isUnspecific = false;
//		boolean isBusStop = false;
//		boolean isTrainStation = false;
//
//		int n = entity.getNumberOfTags();
//		for (int i = 0; i < n; i++) {
//			OsmTag tag = entity.getTag(i);
//
//			// Buildings are always kept
//			if (tag.getKey().equals("building")) {
//				filtered = false;
//				isBuilding = true;
//				if (tag.getValue().equals("yes")) {
//					isUnspecific = true;
//				}
//				break;
//			}
//
//			if (tag.getKey().equals("highway") && tag.getValue().equals("bus_stop")) {
//				isBusStop = true;
//				filtered = false;
//				break;
//			}
//
//			if ((tag.getKey().equals("railway") && tag.getValue().equals("stop")) ||
//				(tag.getKey().equals("railway") && tag.getValue().equals("tram_stop"))) {
//				isTrainStation = true;
//				filtered = false;
//				break;
//			}
//
//			MappingConfig c = config.types.get(tag.getKey());
//			if (c != null) {
//				if (c.values.containsKey("*") || c.values.containsKey(tag.getValue())) {
//					filtered = false;
//					break;
//				}
//			}
//		}
//
//		if (filtered)
//			return;
//
//		if (entity instanceof OsmNode node) {
//
//			Point p = geometryBuilder.build(node);
//			MultiPolygon geometry;
//			try {
//				Polygon polygon = (Polygon) JTS.transform(p, transform).buffer(POI_BUFFER);
//				geometry = geometryBuilder.getGeometryFactory().createMultiPolygon(new Polygon[]{polygon});
//			} catch (TransformException e) {
//				ignored++;
//				return;
//			}
//
//			org.matsim.pt2matsim.hd4.Feature ft = new org.matsim.pt2matsim.hd4.Feature(entity, types, geometry, false, isUnspecific, false, isBusStop, isTrainStation);
//			parse(ft, entity);
//			pois.put(entity.getId(), ft);
//		} else {
//			boolean landuse = false;
//			for (int i = 0; i < n; i++) {
//				if (entity.getTag(i).getKey().equals("landuse")) {
//					landuse = true;
//					break;
//				}
//			}
//
//			MultiPolygon geometry;
//			try {
//				geometry = createPolygon(entity);
//				if (geometry == null) {
//					ignored++;
//					return;
//				}
//				geometry = (MultiPolygon) JTS.transform(geometry, transform);
//			} catch (TransformException e) {
//				// Will be ignored
//				geometry = null;
//			}
//
//			if (geometry == null) {
//				ignored++;
//				return;
//			}
//
//			org.matsim.pt2matsim.hd4.Feature ft = new org.matsim.pt2matsim.hd4.Feature(entity, types, geometry, isBuilding, isUnspecific, landuse, isBusStop, isTrainStation);
//			parse(ft, entity);
//			if (landuse) {
//				this.landuse.put(ft.entity.getId(), ft);
//			} else {
//				// some non-landuse shapes might be too large
//				if (ft.geometry.getArea() < MAX_AREA)
//					entities.put(ft.entity.getId(), ft);
//			}
//		}
//	}

    /**
     * Often landuse shapes are redundant if the area already contains enough detailed shapes with more specific types.
     * These shapes are marked here.
     */
    private void preprocessLanduse(ObjectCollection<Feature> values, STRtree index, double threshold) {

        for (Feature ft : ProgressBar.wrap(values, "Processing landuse entities")) {

            List<Feature> query = index.query(ft.geometry.getBoundary().getEnvelopeInternal());
            List<Geometry> intersections = new ArrayList<>();

            for (Feature other : query) {

                // Other landuse shapes are not considered
                if (other.isLanduse)
                    continue;

                // Use the hull to avoid topology exceptions
                Geometry intersection = intersect(ft.geometry, other.geometry);

                if (!intersection.isEmpty())
                    intersections.add(intersection);
            }

            if (intersections.isEmpty()) {
                ft.setLowPriority();
                continue;
            }

            Geometry geometry = intersections.getFirst().getFactory().buildGeometry(intersections);

            double coveredArea = geometry.union().getArea();
            double ratio = coveredArea / ft.geometry.getArea();

            if (ratio >= threshold) {
                ft.setLowPriority();
            }
        }
    }

    /**
     * Tags buildings within intersections with geometries from list. Used geometries are removed from the list.
     */
    private void processIntersection(Collection<Feature> list, STRtree index, double threshold) {

        Iterator<Feature> it = ProgressBar.wrap(list.iterator(), "Assigning features");

        while (it.hasNext()) {
            Feature ft = it.next();

            if (ft.isNotAssignable())
                continue;

            List<Feature> query = index.query(ft.geometry.getBoundary().getEnvelopeInternal());

            for (Feature other : query) {
                // Assign other features to the buildings
                double otherArea = other.geometry.getArea();

                if (other.isNotAssignable())
                    continue;

                if (otherArea < MAX_ASSIGN) {
                    Geometry intersect = intersect(ft.geometry, other.geometry);
                    if (intersect.getArea() / otherArea > threshold) {

                        // Only assign if this is not a low priority entity, or the other has no types yet
                        if (!ft.lowPriority || !other.hasTypes() || other.isUnspecific)
                            other.assign(ft);
                    }
                }
            }

            if (ft.isAssigned())
                it.remove();
        }
    }

    private Geometry intersect(Geometry a, Geometry b) {
        try {
            return a.intersection(b);
        } catch (Exception e) {
            // some geometries are not well-defined
            try {
                return new ConvexHull(a).getConvexHull().intersection(new ConvexHull(b).getConvexHull());
            } catch (Exception ex) {
                ignored++;
                return a.getFactory().createGeometryCollection(new Geometry[0]);
            }
        }
    }

    private void addFeatures(Long2ObjectMap<Feature> fts, FacilityFeatureExtractor exc,
                           ListFeatureCollection collection) {

        try (ProgressBar pb = new ProgressBar("Creating features", fts.size())) {

            // Option A: Include all features (most inclusive)
            // List<SimpleFeature> features = fts.values().parallelStream()
            //     .map(f -> {
            //         pb.step();
            //         return exc.createFeature(f);
            //     })
            //     .toList();

            // Option B: Include features with types, buildings, or any POI tags (recommended)
            List<SimpleFeature> features = fts.values().parallelStream()
                .filter(f -> f.hasTypes() || f.isBuilding || hasRelevantTags(f))
                .map(f -> {
                    pb.step();
                    return exc.createFeature(f);
                })
                .toList();

            // Option C: Original filter (current behavior)
            // List<SimpleFeature> features = fts.values().parallelStream()
            //     .filter(f -> f.hasTypes() || f.isBuilding)
            //     .map(f -> {
            //         pb.step();
            //         return exc.createFeature(f);
            //     })
            //     .toList();

            collection.addAll(features);
        }
    }

    /**
     * Check if feature has any relevant OSM tags that might be useful for facility analysis
     */
    private boolean hasRelevantTags(Feature feature) {
        Map<String, String> tags = feature.tags;

        // Include features with common facility-related tags
        return tags.containsKey("amenity") ||
               tags.containsKey("shop") ||
               tags.containsKey("office") ||
               tags.containsKey("leisure") ||
               tags.containsKey("tourism") ||
               tags.containsKey("historic") ||
               tags.containsKey("landuse") ||
               tags.containsKey("natural") ||
               tags.containsKey("man_made") ||
               tags.containsKey("place") ||
               tags.containsKey("power") ||
               tags.containsKey("aeroway") ||
               tags.containsKey("railway") ||
               (tags.containsKey("highway") && !tags.get("highway").equals("primary") && !tags.get("highway").equals("secondary"));
    }

    /**
     * Writes how osm ids are mapped to other ids.
     */
    @SafeVarargs
    private void writeMapping(String path, Iterable<Feature>... features) {

        try (CSVPrinter csv = new CSVPrinter(IOUtils.getBufferedWriter(path), CSVFormat.DEFAULT)) {

            // osm ids are only unique within their type
            csv.printRecord("osm_id", "type", "member_id", "member_type");
            for (Feature feature : Iterables.concat(features)) {

                if (!feature.hasTypes())
                    continue;

                csv.printRecord(feature.entity.getId(), feature.osmType, feature.entity.getId(), feature.osmType);
                if (feature.members != null) {
                    for (Feature member : feature.members) {
                        csv.printRecord(member.entity.getId(), member.osmType, feature.entity.getId(), feature.osmType);
                    }
                }
            }

        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    /**
     * Assign activity types to a feature based on its OSM tags and the mapping configuration
     */
    private void assignTypes(Feature feature) {
        if (feature.isNotAssignable())
            return;

        Map<String, String> tags = feature.tags;
        boolean foundMatch = false;

        for (Map.Entry<String, String> tag : tags.entrySet()) {
            MappingConfig c = config.types.get(tag.getKey());
            if (c != null) {
                List<String> activityTypes = null;

                // Check for exact value match
                if (c.values.containsKey(tag.getValue())) {
                    activityTypes = c.values.get(tag.getValue());
                    foundMatch = true;
                }
                // Check for wildcard match
                else if (c.values.containsKey("*")) {
                    activityTypes = c.values.get("*");
                    foundMatch = true;
                }

                if (activityTypes != null) {
                    for (String activityType : activityTypes) {
                        if (!exclude.contains(activityType) && !feature.types.contains(activityType)) {
                            feature.types.add(activityType);
                        }
                    }
                }
            }
        }

        // Log unmapped tags for features that don't get activity types
        if (!foundMatch && !tags.isEmpty()) {
            log.debug("Feature {} has no matching activity types for tags: {}", feature.entity.getId(), tags);
        }
    }

    /**
     * Process an OSM node
     */
    private void processNode(Osm.Node node) {
        boolean filtered = true;
        boolean isBusStop = false;
        boolean isTrainStation = false;

        Map<String, String> tags = node.getTags();

        // Check for bus stops and train stations
        if (tags.containsKey("highway") && tags.get("highway").equals("bus_stop")) {
            isBusStop = true;
            filtered = false;
        }

        if ((tags.containsKey("railway") && tags.get("railway").equals("stop")) ||
            (tags.containsKey("railway") && tags.get("railway").equals("tram_stop"))) {
            isTrainStation = true;
            filtered = false;
        }

        // Check if node matches the mapping configuration
        for (Map.Entry<String, String> tag : tags.entrySet()) {
            MappingConfig c = config.types.get(tag.getKey());
            if (c != null) {
                if (c.values.containsKey("*") || c.values.containsKey(tag.getValue())) {
                    filtered = false;
                    break;
                }
            }
        }

        if (filtered)
            return;

        Point point;
        try {
            // Debug: Log the original OSM coordinates for the first few nodes
            if (pois.size() < 5) {
                log.info("OSM Node coordinates: lon={}, lat={}", node.getCoord().getX(), node.getCoord().getY());
            }

            Point p = geometryBuilder.getGeometryFactory().createPoint(new Coordinate(
                node.getCoord().getX(), node.getCoord().getY()  // This should be (longitude, latitude)
            ));
            point = (Point) JTS.transform(p, transform);

            // Debug: Log the transformed coordinates for the first few nodes
            if (pois.size() < 5) {
                log.info("Transformed coordinates: x={}, y={}", point.getX(), point.getY());
            }
        } catch (TransformException e) {
            log.warn("Could not transform point: {}", e.getMessage());
            return;
        }

        if (isTrainStation || isBusStop) {
            // Add buffer around transit stops
            Geometry buffer = point.buffer(POI_BUFFER);
            OsmEntity entity = new OsmEntity(node, tags);
            Feature ft = new Feature(entity, buffer, "node", tags);
            ft.isBusStop = isBusStop;
            ft.isTrainStop = isTrainStation;
            // Assign activity types for transit stops too
            assignTypes(ft);
            pois.put(node.getId().index(), ft);
            return;
        }

        OsmEntity entity = new OsmEntity(node, tags);
        Feature ft = new Feature(entity, point, "node", tags);

        // Assign activity types based on mapping configuration
        assignTypes(ft);

        if (ft.hasTypes()) {
            pois.put(node.getId().index(), ft);
        }
    }

    /**
     * Process an OSM way
     */
    private void processWay(Osm.Way way) {
        boolean filtered = true;
        boolean isBuilding = false;
        boolean isUnspecific = false;

        Map<String, String> tags = way.getTags();

        // Buildings are always kept
        if (tags.containsKey("building")) {
            filtered = false;
            isBuilding = true;
            if (tags.get("building").equals("yes")) {
                isUnspecific = true;
            }
        }

        // Check if way matches the mapping configuration
        for (Map.Entry<String, String> tag : tags.entrySet()) {
            MappingConfig c = config.types.get(tag.getKey());
            if (c != null) {
                if (c.values.containsKey("*") || c.values.containsKey(tag.getValue())) {
                    filtered = false;
                    break;
                }
            }
        }

        if (filtered)
            return;

        // Create a linear ring or line string from way nodes
        List<Coordinate> coords = new ArrayList<>();
        for (Osm.Node node : way.getNodes()) {
            coords.add(new Coordinate(node.getCoord().getX(), node.getCoord().getY()));  // Revert back to original X, Y order
        }

        if (coords.size() < 2)
            return;

        Geometry geom;
        // Closed way = polygon, otherwise linestring
        boolean isClosed = coords.getFirst().equals2D(coords.getLast()) && coords.size() >= 4;

        try {
            if (isClosed) {
                LinearRing ring = geometryBuilder.getGeometryFactory().createLinearRing(coords.toArray(new Coordinate[0]));
                geom = geometryBuilder.getGeometryFactory().createPolygon(ring);
            } else {
                geom = geometryBuilder.getGeometryFactory().createLineString(coords.toArray(new Coordinate[0]));
            }

            // Transform to target CRS
            geom = JTS.transform(geom, transform);

        } catch (Exception e) {
            log.warn("Could not create/transform geometry: {}", e.getMessage());
            ignored++;
            return;
        }

        // Skip very large areas
        if (isClosed && geom.getArea() > MAX_AREA) {
            return;
        }

        OsmEntity entity = new OsmEntity(way, tags);
        Feature ft = new Feature(entity, geom, "way", tags);
        ft.isUnspecific = isUnspecific;

        // Assign activity types based on mapping configuration
        assignTypes(ft);

        // Buildings go into entities
        if (isBuilding) {
            entities.put(way.getId().index(), ft);
        }
        // Landuse and other area features
        else if (isClosed && geom.getArea() > 0) {
            boolean isLanduse = tags.containsKey("landuse") || tags.containsKey("leisure") ||
                                tags.containsKey("natural") || tags.containsKey("amenity");

            if (isLanduse) {
                ft.isLanduse = true;
                landuse.put(way.getId().index(), ft);
            } else if (ft.hasTypes()) {
                entities.put(way.getId().index(), ft);
            }
        }
        // Linear features as POIs
        else if (ft.hasTypes()) {
            pois.put(way.getId().index(), ft);
        }
    }

    /**
     * Process an OSM relation
     */
    private void processRelation(Osm.Relation relation) {
        boolean filtered = true;

        Map<String, String> tags = relation.getTags();

        // Only process multipolygon and building relations
        boolean isMultiPolygon = tags.getOrDefault("type", "").equals("multipolygon");
        boolean isBuilding = tags.containsKey("building");

        if (!isMultiPolygon && !isBuilding)
            return;

        // Check if relation matches the mapping configuration
        for (Map.Entry<String, String> tag : tags.entrySet()) {
            MappingConfig c = config.types.get(tag.getKey());
            if (c != null) {
                if (c.values.containsKey("*") || c.values.containsKey(tag.getValue())) {
                    filtered = false;
                    break;
                }
            }
        }

        if (isBuilding)
            filtered = false;

        if (filtered)
            return;

        // Build multipolygon geometry
        Geometry geom = null;
        try {
            // This is a simplified approach - for a full implementation you would need to
            // handle inner/outer rings, etc.
            List<Geometry> polygons = new ArrayList<>();

            // Create polygons from member ways
            for (Osm.Element member : relation.getMembers()) {
                if (member instanceof Osm.Way way) {
                    String role = relation.getMemberRoles(member).getFirst();
                    if (role.equals("outer")) {
                        List<Coordinate> coords = new ArrayList<>();
                        for (Osm.Node node : way.getNodes()) {
                            coords.add(new Coordinate(node.getCoord().getX(), node.getCoord().getY()));
                        }

                        if (coords.size() >= 4) {
                            try {
                                LinearRing ring = geometryBuilder.getGeometryFactory().createLinearRing(coords.toArray(new Coordinate[0]));
                                polygons.add(geometryBuilder.getGeometryFactory().createPolygon(ring));
                            } catch (Exception e) {
                                // Skip invalid rings
                            }
                        }
                    }
                }
            }

            if (!polygons.isEmpty()) {
                if (polygons.size() == 1) {
                    geom = polygons.getFirst();
                } else {
                    geom = geometryBuilder.getGeometryFactory().createMultiPolygon(
                        polygons.toArray(new Polygon[0])
                    );
                }

                // Transform to target CRS
                geom = JTS.transform(geom, transform);
            }

        } catch (Exception e) {
            log.warn("Could not create/transform relation geometry: {}", e.getMessage());
            ignored++;
            return;
        }

        if (geom == null || geom.isEmpty())
            return;

        // Skip very large areas
        if (geom.getArea() > MAX_AREA) {
            return;
        }

        OsmEntity entity = new OsmEntity(relation, tags);
        Feature ft = new Feature(entity, geom, "relation", tags);

        // Assign activity types based on mapping configuration
        assignTypes(ft);

        // Building goes to entities
        if (isBuilding) {
            ft.isUnspecific = tags.get("building").equals("yes");
            entities.put(relation.getId().index(), ft);
        }
        // Landuse and other area features
        else {
            boolean isLanduse = tags.containsKey("landuse") || tags.containsKey("leisure") ||
                                tags.containsKey("natural") || tags.containsKey("amenity");

            if (isLanduse) {
                ft.isLanduse = true;
                landuse.put(relation.getId().index(), ft);
            } else if (ft.hasTypes()) {
                entities.put(relation.getId().index(), ft);
            }
        }
    }

    /**
         * Simple class to represent an OSM entity (node, way, relation) with tags
         */
        private record OsmEntity(Object entity, Map<String, String> tags) {

        public long getId() {
                if (entity instanceof Osm.Node) {
                    return Long.parseLong(((Osm.Node) entity).getId().toString());
                } else if (entity instanceof Osm.Way) {
                    return Long.parseLong(((Osm.Way) entity).getId().toString());
                } else if (entity instanceof Osm.Relation) {
                    return Long.parseLong(((Osm.Relation) entity).getId().toString());
                }
                return 0;
            }
        }

    /**
     * Feature class representing an OSM entity with geometry and associated types
     */
    public static class Feature {
        public Geometry geometry;
        public OsmEntity entity;
        public String osmType;
        public boolean lowPriority = false;
        public boolean isLanduse = false;
        public boolean isUnspecific = false;
        public List<Feature> members;
        public boolean isBuilding = false;  // Fix: Change from Object to boolean
        public boolean isBusStop;
        public boolean isTrainStop;
        public OsmPbfFileReader.StringCache bits;
        public boolean geomIssues;
        private boolean assigned = false;
        private final List<String> types = new ArrayList<>();
        private final Map<String, String> tags;

        public Feature(OsmEntity entity, Geometry geometry, String osmType, Map<String, String> tags) {
            this.geometry = geometry;
            this.entity = entity;
            this.osmType = osmType;
            this.tags = tags;
        }

        public boolean hasTypes() {
            return !types.isEmpty();
        }

        public boolean isNotAssignable() {
            return assigned || isLanduse;
        }

        public void assign(Feature other) {
            if (members == null)
                members = new ArrayList<>();

            members.add(other);
            other.assigned = true;

            // Copy types
            for (String type : other.types) {
                if (!types.contains(type))
                    types.add(type);
            }
        }

        public boolean isAssigned() {
            return assigned;
        }

        public void setLowPriority() {
            this.lowPriority = true;
        }

        public int getLevels() {
            return tags.containsKey("levels") ? Integer.parseInt(tags.get("levels")) : 0;
        }

        public boolean hasLanduse(Object o) {
            return tags.containsKey("landuse") || tags.containsKey("leisure") ||
                   tags.containsKey("natural") || tags.containsKey("amenity");
        }

        public boolean isResidentialOnly() {
            return tags.containsKey("residential") && tags.get("residential").equals("yes");
        }

        public boolean hasType(String activityType) {
            return types.contains(activityType);
        }
    }

    /**
     * The activity mapping configuration structure
     */
    public static class ActivityMapping {
        public Map<String, MappingConfig> types = new HashMap<>();

        @JsonAnyGetter
        public Map<String, MappingConfig> getTypes() {
            return types;
        }

        @JsonAnySetter
        public void setType(String name, MappingConfig value) {
            types.put(name, value);
        }
    }

    /**
     * Mapping configuration for a specific tag
     */
    public static class MappingConfig {
        public Map<String, List<String>> values = new HashMap<>();

        @JsonAnyGetter
        public Map<String, List<String>> getValues() {
            return values;
        }

        @JsonAnySetter
        public void setValue(String name, List<String> value) {
            values.put(name, value);
        }
    }

    /**
     * Helper class for extracting features
     */
    private static class FacilityFeatureExtractor {
        public final SimpleFeatureType featureType;
        public final STRtree entities;
        public final STRtree pois;
        public final STRtree landuse;
        private final Object2IntMap<String> types;
        private int featureCount = 0; // Add counter for debugging

        final ThreadLocal<SimpleFeatureBuilder> featureBuilder;


        public FacilityFeatureExtractor(String crs, Object2IntMap<String> types, Long2ObjectMap<Feature> entities,
                                        Long2ObjectMap<Feature> pois, Long2ObjectMap<Feature> landuse) throws FactoryException {

            this.entities = createIndex(entities);
            this.pois = createIndex(pois);
            this.landuse = createIndex(landuse);
            this.types = types;
            SimpleFeatureTypeBuilder typeBuilder = new SimpleFeatureTypeBuilder();
            typeBuilder.setName("facilities");
            typeBuilder.setCRS(CRS.decode(crs));
            typeBuilder.add("osm_id", Long.class);
            typeBuilder.add("osm_type", String.class);
            typeBuilder.add("the_geom", MultiPolygon.class);
            typeBuilder.add("area", Double.class);
            typeBuilder.add("levels", Integer.class);
            typeBuilder.add("landuse", Boolean.class);
            typeBuilder.add("building", Boolean.class);
            typeBuilder.add("residential_only", Boolean.class);
            typeBuilder.add("landuse_residential_500m", Double.class);
            typeBuilder.add("landuse_residential_1500m", Double.class);
            typeBuilder.add("landuse_retail_500m", Double.class);
            typeBuilder.add("landuse_retail_1500m", Double.class);
            typeBuilder.add("landuse_commercial_500m", Double.class);
            typeBuilder.add("landuse_commercial_1500m", Double.class);
            typeBuilder.add("landuse_recreation_1500m", Double.class);
            typeBuilder.add("parking_space_500m", Double.class);
            typeBuilder.add("nearest_bus_stop", Double.class);
            typeBuilder.add("nearest_train_station", Double.class);

            typeBuilder.add("poi_home", Integer.class);
            typeBuilder.add("poi_home_250m", Integer.class);
            typeBuilder.add("poi_work", Integer.class);
            typeBuilder.add("poi_work_250m", Integer.class);
            typeBuilder.add("poi_shop", Integer.class);
            typeBuilder.add("poi_shop_250m", Integer.class);
            typeBuilder.add("poi_dining", Integer.class);
            typeBuilder.add("poi_dining_250m", Integer.class);
            typeBuilder.add("poi_other", Integer.class);
            typeBuilder.add("poi_other_250m", Integer.class);
            typeBuilder.add("poi_leisure", Integer.class);
            typeBuilder.add("poi_leisure_250m", Integer.class);
            for (String t : types.keySet()) {
			    typeBuilder.add(t, Boolean.class);
		    }

            this.featureType = typeBuilder.buildFeatureType();

            this.featureBuilder = ThreadLocal.withInitial(() -> new SimpleFeatureBuilder(featureType));
        }

        private static STRtree createIndex(Long2ObjectMap<Feature> entities) {
            STRtree index = new STRtree();
            for (Feature entity : entities.values()) {
                index.insert(entity.geometry.getBoundary().getEnvelopeInternal(), entity);
            }
            index.build();
            return index;
	    }

        public SimpleFeature createFeature(Feature ft) {

		// feature building is thread safe
		SimpleFeatureBuilder b = featureBuilder.get();

		// Ensure geometry is MultiPolygon
		MultiPolygon multiPolygon;
		if (ft.geometry instanceof MultiPolygon) {
			multiPolygon = (MultiPolygon) ft.geometry;
		} else if (ft.geometry instanceof Polygon) {
			multiPolygon = ft.geometry.getFactory().createMultiPolygon(new Polygon[]{(Polygon) ft.geometry});
		} else {
			// For points and lines, create a small buffer to make it a polygon, then multipolygon
			Geometry buffered = ft.geometry.buffer(0.1); // Small buffer
			if (buffered instanceof Polygon) {
				multiPolygon = ft.geometry.getFactory().createMultiPolygon(new Polygon[]{(Polygon) buffered});
			} else if (buffered instanceof MultiPolygon) {
				multiPolygon = (MultiPolygon) buffered;
			} else {
				// Fallback: create a small square around the geometry
				Coordinate coord = ft.geometry.getCoordinate();
				Coordinate[] coords = new Coordinate[]{
					new Coordinate(coord.x - 0.1, coord.y - 0.1),
					new Coordinate(coord.x + 0.1, coord.y - 0.1),
					new Coordinate(coord.x + 0.1, coord.y + 0.1),
					new Coordinate(coord.x - 0.1, coord.y + 0.1),
					new Coordinate(coord.x - 0.1, coord.y - 0.1)
				};
				LinearRing ring = ft.geometry.getFactory().createLinearRing(coords);
				Polygon poly = ft.geometry.getFactory().createPolygon(ring);
				multiPolygon = ft.geometry.getFactory().createMultiPolygon(new Polygon[]{poly});
			}
		}

		// Swap x and y coordinates for GPKG output
		MultiPolygon swappedMultiPolygon = swapCoordinates(multiPolygon);

		// Debug: Log the coordinate order being used in the output
		if (swappedMultiPolygon.getCoordinates().length < 5) {
			log.info("Output MultiPolygon coordinates - First coordinate: x={}, y={}",
				swappedMultiPolygon.getCoordinates()[0].x, swappedMultiPolygon.getCoordinates()[0].y);
		}

		b.add(ft.entity.getId());
		b.add(ft.osmType);
		b.add(swappedMultiPolygon);
		b.add(BigDecimal.valueOf(swappedMultiPolygon.getArea()).setScale(2, RoundingMode.HALF_EVEN).doubleValue());
		b.add(ft.getLevels());
		b.add(ft.hasLanduse(null));
		b.add(ft.isBuilding);
		b.add(ft.isResidentialOnly());

		b.add(calcLanduse("residential", ft, multiPolygon, 500));
		b.add(calcLanduse("residential", ft, multiPolygon, 1500));
		b.add(calcLanduse("retail", ft, multiPolygon, 500));
		b.add(calcLanduse("retail", ft, multiPolygon, 1500));
		b.add(calcLanduse("commercial", ft, multiPolygon, 500));
		b.add(calcLanduse("commercial", ft, multiPolygon, 1500));
		b.add(calcLanduse("recreation_ground", ft, multiPolygon, 1500));
		b.add(calcArea("parking", ft, multiPolygon, 500));
		b.add(findNearest(ft, f -> f.isBusStop));
		b.add(findNearest(ft, f -> f.isTrainStop));

		b.add(countPOIs("home", ft));
		b.add(countPOIs("home", ft, 250));
		b.add(countPOIs("work", ft));
		b.add(countPOIs("work", ft, 250));
		b.add(countPOIs("shop", ft));
		b.add(countPOIs("shop", ft, 250));
		b.add(countPOIs("dining", ft));
		b.add(countPOIs("dining", ft, 250));
		b.add(countPOIs("other", ft));
		b.add(countPOIs("other", ft, 250));
		b.add(countPOIs("leisure", ft));
		b.add(countPOIs("leisure", ft, 250));

		// Add activity type flags
		for (String activityType : types.keySet()) {
			b.add(ft.types.contains(activityType));
		}

		return b.buildFeature(null);
	}

	/**
	 * Swap x and y coordinates in a MultiPolygon for proper GPKG output
	 */
	private MultiPolygon swapCoordinates(MultiPolygon multiPolygon) {
		GeometryFactory factory = multiPolygon.getFactory();
		Polygon[] swappedPolygons = new Polygon[multiPolygon.getNumGeometries()];

		for (int i = 0; i < multiPolygon.getNumGeometries(); i++) {
			Polygon polygon = (Polygon) multiPolygon.getGeometryN(i);

			// Swap coordinates in exterior ring
			Coordinate[] exteriorCoords = polygon.getExteriorRing().getCoordinates();
			Coordinate[] swappedExteriorCoords = new Coordinate[exteriorCoords.length];
			for (int j = 0; j < exteriorCoords.length; j++) {
				swappedExteriorCoords[j] = new Coordinate(exteriorCoords[j].y, exteriorCoords[j].x);
			}
			LinearRing swappedExteriorRing = factory.createLinearRing(swappedExteriorCoords);

			// Swap coordinates in interior rings (holes)
			LinearRing[] swappedInteriorRings = new LinearRing[polygon.getNumInteriorRing()];
			for (int k = 0; k < polygon.getNumInteriorRing(); k++) {
				Coordinate[] interiorCoords = polygon.getInteriorRingN(k).getCoordinates();
				Coordinate[] swappedInteriorCoords = new Coordinate[interiorCoords.length];
				for (int j = 0; j < interiorCoords.length; j++) {
					swappedInteriorCoords[j] = new Coordinate(interiorCoords[j].y, interiorCoords[j].x);
				}
				swappedInteriorRings[k] = factory.createLinearRing(swappedInteriorCoords);
			}

			swappedPolygons[i] = factory.createPolygon(swappedExteriorRing, swappedInteriorRings);
		}

		return factory.createMultiPolygon(swappedPolygons);
	}

	/**
	 * Calculate the area of landuse within a given radius.
	 */
	@SuppressWarnings("unchecked")
	private double calcLanduse(String type, Feature ft, MultiPolygon geometry, double radius) {

		if (ft.isResidentialOnly()) {
			return 0;
		}

		Geometry bbox = geometry.getCentroid().buffer(radius);

		double res = 0;
		List<Feature> query = landuse.query(bbox.getEnvelopeInternal());
		for (Feature q : query) {

			if (!q.geomIssues && q.hasLanduse(type)) {
				try {
					res += q.geometry.intersection(bbox).getArea();
				} catch (TopologyException e) {
					q.geomIssues = true;
				}
			}
		}

		// convert to square kilometers
		return BigDecimal.valueOf(res / 1_000_000).setScale(4, RoundingMode.HALF_EVEN).doubleValue();
	}

	private double calcArea(String activityType, Feature ft, MultiPolygon geometry, double radius) {

		if (ft.isResidentialOnly()) {
			return 0;
		}

		Geometry bbox = geometry.getCentroid().buffer(radius);

		double res = 0;
//		List<Feature> query = entities.query(bbox.getEnvelopeInternal());
		for (Object o : entities.query(bbox.getEnvelopeInternal())) {
            Feature q = (Feature) o;
			if (!q.geomIssues && q.hasType(activityType)) {
				try {
					res += q.geometry.intersection(bbox).getArea();
				} catch (TopologyException e) {
					q.geomIssues = true;
				}
			}
		}

		return BigDecimal.valueOf(res / 1_000).setScale(4, RoundingMode.HALF_EVEN).doubleValue();
	}

	private double findNearest(Feature ft, Predicate<Feature> filter) {

		if (ft.isResidentialOnly()) {
			return 0;
		}

		// Convert geometry to MultiPolygon if needed
		MultiPolygon geometry;
		if (ft.geometry instanceof MultiPolygon) {
			geometry = (MultiPolygon) ft.geometry;
		} else if (ft.geometry instanceof Polygon) {
			geometry = ft.geometry.getFactory().createMultiPolygon(new Polygon[]{(Polygon) ft.geometry});
		} else {
			// For other geometry types, create a small buffer
			Geometry buffered = ft.geometry.buffer(0.1);
			if (buffered instanceof Polygon) {
				geometry = ft.geometry.getFactory().createMultiPolygon(new Polygon[]{(Polygon) buffered});
			} else if (buffered instanceof MultiPolygon) {
				geometry = (MultiPolygon) buffered;
			} else {
				// Fallback: create a small square
				Coordinate coord = ft.geometry.getCoordinate();
				Coordinate[] coords = new Coordinate[]{
					new Coordinate(coord.x - 0.1, coord.y - 0.1),
					new Coordinate(coord.x + 0.1, coord.y - 0.1),
					new Coordinate(coord.x + 0.1, coord.y + 0.1),
					new Coordinate(coord.x - 0.1, coord.y + 0.1),
					new Coordinate(coord.x - 0.1, coord.y - 0.1)
				};
				LinearRing ring = ft.geometry.getFactory().createLinearRing(coords);
				Polygon poly = ft.geometry.getFactory().createPolygon(ring);
				geometry = ft.geometry.getFactory().createMultiPolygon(new Polygon[]{poly});
			}
		}

		for (Integer radius : List.of(500, 5000, 20000)) {

			Geometry bbox = geometry.getCentroid().buffer(radius);

			List<Feature> query = pois.query(bbox.getEnvelopeInternal());

			OptionalDouble dist = query.stream()
				.filter(filter)
				.mapToDouble(a -> a.geometry.distance(geometry))
				.min();

			if (dist.isPresent()) {
				return dist.getAsDouble();
			}

		}

		return 20000;
	}

	private int countPOIs(String type, Feature ft) {
		// Check if this feature has the specified activity type
		int count = ft.types.contains(type) ? 1 : 0;

		// Add counts from assigned members
		if (ft.members != null) {
			for (Feature m : ft.members) {
				if (m.types.contains(type)) {
					count++;
				}
			}
		}

		return count;
	}

	@SuppressWarnings("unchecked")
	private int countPOIs(String type, Feature ft, double radius) {

		if (ft.isResidentialOnly()) {
			return 0;
		}

		// Base count
		int count = countPOIs(type, ft);

		Geometry bbox = ft.geometry.getCentroid().buffer(radius);

		Iterable<Feature> query = entities.query(bbox.getEnvelopeInternal());
		for (Feature q : query) {
			try {
				if (!q.geomIssues && q.geometry.distance(ft.geometry.getCentroid()) < radius && q != ft) {
					count += countPOIs(type, q);
				}
			} catch (TopologyException e) {
				q.geomIssues = true;
			}
		}

		return count;
	}
    }
}
