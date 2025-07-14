package org.matsim.pt2matsim.hd4;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.geotools.api.data.FileDataStore;
import org.geotools.api.data.FileDataStoreFinder;
import org.geotools.api.data.SimpleFeatureSource;
import org.geotools.api.feature.simple.SimpleFeature;
import org.geotools.api.referencing.FactoryException;
import org.geotools.feature.FeatureIterator;
import org.locationtech.jts.geom.*;
import org.locationtech.jts.index.strtree.STRtree;
import org.matsim.api.core.v01.Coord;
import org.matsim.api.core.v01.network.Link;
import org.matsim.api.core.v01.network.Network;
import org.matsim.application.MATSimAppCommand;
import org.matsim.application.options.ShpOptions;
import org.matsim.core.utils.io.IOUtils;
import org.matsim.pt.transitSchedule.api.TransitSchedule;
import org.matsim.pt2matsim.config.OsmConverterConfigGroup;
import org.matsim.pt2matsim.config.PublicTransitMappingConfigGroup;
import org.matsim.pt2matsim.gtfs.GtfsFeed;
import org.matsim.pt2matsim.mapping.PTMapper;
import org.matsim.pt2matsim.run.CheckMappedSchedulePlausibility;
import org.matsim.pt2matsim.run.Osm2MultimodalNetwork;
import org.matsim.pt2matsim.run.gis.Network2Geojson;
import org.matsim.pt2matsim.run.gis.Schedule2Geojson;
import org.matsim.pt2matsim.tools.NetworkTools;
import org.matsim.pt2matsim.tools.ScheduleTools;
import picocli.CommandLine;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

/**
 * Workflow_SG class handles the workflow for processing and mapping public transit data in Singapore.
 */
@CommandLine.Command(name = "Workflow_SG", description = "Workflow for processing and mapping public transit data in Singapore.")
public class Workflow_SG_HD4 implements MATSimAppCommand {
    private static final Logger log = LogManager.getLogger(Workflow_SG_HD4.class);

    // Directory paths and file names
    @CommandLine.Option(names = "--base-dir", description = "Base directory path", defaultValue = "")
    public static String BASE_DIR;
    
    @CommandLine.Option(names = "--output-dir", description = "Output directory path", defaultValue = "output/islandwide/")
    private static String OUTPUT_DIR;
    
    @CommandLine.Option(names = "--gtfs-dir", description = "GTFS directory path", defaultValue = "resources/GTFS/gtfs_2022/")
    private String GTFS_DIR;
    
    @CommandLine.Option(names = "--unmapped-schedule-bus", description = "Unmapped bus schedule file path", defaultValue = "sg-2023_unmapped_transitSchedule_epsg3414_bus.xml.gz")
    private String UNMAPPED_SCHEDULE_BUS;
    
    @CommandLine.Option(names = "--unmapped-schedule-train", description = "Unmapped train schedule file path", defaultValue ="resources/sg-2023_unmapped_transitSchedule_epsg3414_train.xml")
    private String UNMAPPED_SCHEDULE_TRAIN;
    
    @CommandLine.Option(names = "--unmapped-schedule-merge", description = "Unmapped merged schedule file path", defaultValue = "sg-2023_unmapped_transitSchedule_epsg3414_merge.xml.gz")
    private static String UNMAPPED_SCHEDULE_MERGE;
    
    @CommandLine.Option(names = "--temp-network-xml", description = "Temporary network XML file path", defaultValue = "sg-2023_network1_epsg3414_temp.xml.gz")
    private static String TEMP_NETWORK_XML;
    
    @CommandLine.Option(names = "--input-osm-file", description = "Input OSM file path", defaultValue ="input/singapore_all_240101.osm")
    private String INPUT_OSM_FILE;
    
    @CommandLine.Option(names = "--final-normal-network-xml", description = "Final normal network XML file path", defaultValue = "sg-2023_network_epsg3414.xml.gz")
    public static String FINAL_NORMAL_NETWORK_XML;
    
    @CommandLine.Option(names = "--final-normal-schedule-xml", description = "Final normal schedule XML file path", defaultValue = "sg-2023_transitSchedule_epsg3414.xml.gz")
    public static String FINAL_NORMAL_SCHEDULE_XML;
    
    @CommandLine.Option(names = "--final-shapes-network-xml", description = "Final shapes network XML file path", defaultValue = "sg-2023_network_epsg3414_shapes.xml.gz")
    private static String FINAL_SHAPES_NETWORK_XML;
    
    @CommandLine.Option(names = "--final-shapes-schedule-xml", description = "Final shapes schedule XML file path", defaultValue = "sg-2023_transitSchedule_epsg3414_shapes.xml.gz")
    private static String FINAL_SHAPES_SCHEDULE_XML;
    
    @CommandLine.Option(names = "--final-additional-line-info-file", description = "Final additional line info file path", defaultValue = "sg-2023_additionalLineInfo.xml.gz")
    private String FINAL_ADDITIONAL_LINE_INFO_FILE;
    
    @CommandLine.Option(names = "--final-vehicle-file-xml", description = "Final vehicle file XML path", defaultValue = "sg-2023_transitVehicles.xml.gz")
    private String FINAL_VEHICLE_FILE_XML;
    
    @CommandLine.Option(names = "--boundary-shp-file", description = "Boundary shapefile path", defaultValue ="input/clementi_PA.shp")
    private String BOUNDARY_SHP_FILE;
    
    @CommandLine.Option(names = "--final-normal-network-geojson", description = "Final normal network GeoJSON file path", defaultValue = "sg-2023_network_epsg3414.geojson")
    private static String FINAL_NORMAL_NETWORK_GEOJSON;
    
    @CommandLine.Option(names = "--final-normal-schedule-geojson", description = "Final normal schedule GeoJSON file path", defaultValue = "sg-2023_transitSchedule_epsg3414.geojson")
    private static String FINAL_NORMAL_SCHEDULE_GEOJSON;
    
    @CommandLine.Option(names = "--final-shapes-network-geojson", description = "Final shapes network GeoJSON file path", defaultValue = "sg-2023_network_epsg3414_shapes.geojson")
    private static String FINAL_SHAPES_NETWORK_GEOJSON;
    
    @CommandLine.Option(names = "--final-shapes-schedule-geojson", description = "Final shapes schedule GeoJSON file path", defaultValue = "sg-2023_transitSchedule_epsg3414_shapes.geojson")
    private static String FINAL_SHAPES_SCHEDULE_GEOJSON;
    
    @CommandLine.Option(names = "--output-detailed-link-geometry-csv-file", description = "Output detailed link geometry file path", defaultValue = "sg-2023_network_epsg3414.csv")
    private String OUTPUT_DETAILED_LINK_GEOMETRY_FILE;
    
    @CommandLine.Option(names = "--output-coordinate-system", description = "Output coordinate system", defaultValue = "EPSG:3414")
    private static String OUTPUT_COORDINATE_SYSTEM;

    private static GtfsFeed gtfsFeed;

    
    @CommandLine.Mixin
    private ShpOptions shpOptions = new ShpOptions();
    
    /**
     * Main method to execute the workflow.
     *
     * @param args Command line arguments.
     * @throws FactoryException If there is an error with the factory.
     * @throws IOException If there is an I/O error.
     */
    public static void main(String[] args) throws FactoryException, IOException {new Workflow_SG_HD4().execute(args);}

    public static String getOutputCRS() {
        return OUTPUT_COORDINATE_SYSTEM;
    }

    @Override
    public Integer call() throws Exception {
        OUTPUT_DIR = BASE_DIR + OUTPUT_DIR;
        GTFS_DIR = BASE_DIR + GTFS_DIR;
        UNMAPPED_SCHEDULE_BUS = OUTPUT_DIR + UNMAPPED_SCHEDULE_BUS;
        UNMAPPED_SCHEDULE_TRAIN = BASE_DIR + UNMAPPED_SCHEDULE_TRAIN;
        UNMAPPED_SCHEDULE_MERGE = OUTPUT_DIR + UNMAPPED_SCHEDULE_MERGE;
        TEMP_NETWORK_XML = OUTPUT_DIR + TEMP_NETWORK_XML;
        INPUT_OSM_FILE = BASE_DIR + INPUT_OSM_FILE;
        FINAL_NORMAL_NETWORK_XML = OUTPUT_DIR + FINAL_NORMAL_NETWORK_XML;
        FINAL_NORMAL_SCHEDULE_XML = OUTPUT_DIR + FINAL_NORMAL_SCHEDULE_XML;
        FINAL_SHAPES_NETWORK_XML = OUTPUT_DIR + FINAL_SHAPES_NETWORK_XML;
        FINAL_SHAPES_SCHEDULE_XML = OUTPUT_DIR + FINAL_SHAPES_SCHEDULE_XML;
        FINAL_ADDITIONAL_LINE_INFO_FILE = OUTPUT_DIR + FINAL_ADDITIONAL_LINE_INFO_FILE;
        FINAL_VEHICLE_FILE_XML = OUTPUT_DIR + FINAL_VEHICLE_FILE_XML;
        BOUNDARY_SHP_FILE = BASE_DIR + BOUNDARY_SHP_FILE;
        FINAL_NORMAL_NETWORK_GEOJSON = OUTPUT_DIR + FINAL_NORMAL_NETWORK_GEOJSON;
        FINAL_NORMAL_SCHEDULE_GEOJSON = OUTPUT_DIR + FINAL_NORMAL_SCHEDULE_GEOJSON;
        FINAL_SHAPES_NETWORK_GEOJSON = OUTPUT_DIR + FINAL_SHAPES_NETWORK_GEOJSON;
        FINAL_SHAPES_SCHEDULE_GEOJSON = OUTPUT_DIR + FINAL_SHAPES_SCHEDULE_GEOJSON;
        OUTPUT_DETAILED_LINK_GEOMETRY_FILE = OUTPUT_DIR + OUTPUT_DETAILED_LINK_GEOMETRY_FILE;
        
        
        // Replace lines 146-148 with this code
        if (!new File(OUTPUT_DIR).exists()) {
            boolean created = new File(OUTPUT_DIR).mkdirs();
            if (created) {
                log.info("Created output directory: {}", OUTPUT_DIR);
            } else {
                log.error("Failed to create output directory: {}", OUTPUT_DIR);
            }
        }

        OsmConverterConfigGroup osmConfig = OsmConverterConfigGroup.loadConfig(BASE_DIR + "input/Osm2MATSimConfig_CARES.xml");
        osmConfig.setOutputNetworkFile(TEMP_NETWORK_XML);
        osmConfig.setOsmFile(INPUT_OSM_FILE);
        osmConfig.setOutputDetailedLinkGeometryFile(OUTPUT_DETAILED_LINK_GEOMETRY_FILE);
        osmConfig.setOutputCoordinateSystem(OUTPUT_COORDINATE_SYSTEM);

        Osm2MultimodalNetwork.run(osmConfig);

        log.info("==============================================================");
        log.info("=      Cutting off network outside boundary of Singapore     =");
        log.info("==============================================================");

        Network network = NetworkTools.readNetwork(osmConfig.getOutputNetworkFile());
        Network2Geojson.run(OUTPUT_COORDINATE_SYSTEM, network, osmConfig.getOutputNetworkFile().replace(".xml.gz", ".geojson"));
        cutOffNetwork(network, BOUNDARY_SHP_FILE);
        NetworkTools.writeNetwork(network, osmConfig.getOutputNetworkFile().replace(".xml.gz", "_cleaned.xml.gz"));
        log.info("==> Network cut off successfully!");
        log.info("==============================================================");
        Network2Geojson.run(OUTPUT_COORDINATE_SYSTEM, network, osmConfig.getOutputNetworkFile().replace(".xml.gz", "_cleaned.geojson"));
//
//        Gtfs2TransitSchedule.run(GTFS_DIR, "dayWithMostTrips", osmConfig.getOutputCoordinateSystem(), UNMAPPED_SCHEDULE_BUS, FINAL_VEHICLE_FILE_XML, FINAL_ADDITIONAL_LINE_INFO_FILE);
//
//        log.info("==============================================================");
//        log.info("=                   Merging with train schedule              =");
//        log.info("==============================================================");
//
//        TransitSchedule schedule = ScheduleTools.readTransitSchedule(UNMAPPED_SCHEDULE_BUS);
//        TransitSchedule scheduleTrain = ScheduleTools.readTransitSchedule(UNMAPPED_SCHEDULE_TRAIN);
//
////        ScheduleTools.removeLinkRefs(scheduleTrain);
//        ScheduleTools.mergeSchedules(schedule, scheduleTrain);
//        ScheduleTools.writeTransitSchedule(schedule, UNMAPPED_SCHEDULE_MERGE);
//
//        log.info("==> Train schedule merged successfully!");
//        log.info("==============================================================");
//
//        mappingAnalysisNormal();
//        mappingAnalysisWithShapes();
        return 0;
    }
    
    /**
     * Creates a configuration for public transit mapping.
     *
     * @return PublicTransitMappingConfigGroup configuration object.
     */
    public static PublicTransitMappingConfigGroup createTestPTMConfig() {
        PublicTransitMappingConfigGroup config = new PublicTransitMappingConfigGroup();

        // Configure transport mode assignments
        PublicTransitMappingConfigGroup.TransportModeAssignment mraBus = new PublicTransitMappingConfigGroup.TransportModeAssignment("bus");
        mraBus.setNetworkModesStr("car,bus");
        config.addParameterSet(mraBus);

        PublicTransitMappingConfigGroup.TransportModeAssignment tmaSubway = new PublicTransitMappingConfigGroup.TransportModeAssignment("subway");
        tmaSubway.setNetworkModesStr("subway");
        config.addParameterSet(tmaSubway);
        PublicTransitMappingConfigGroup.TransportModeAssignment tmaMonorail = new PublicTransitMappingConfigGroup.TransportModeAssignment("monorail");
        tmaMonorail.setNetworkModesStr("monorail");
        config.addParameterSet(tmaMonorail);

        // Set additional configuration parameters
        config.setNumOfThreads(16);
        config.setMaxLinkCandidateDistance(100);
        config.setModesToKeepOnCleanUp(Set.of("car", "bus", "subway", "bike", "walk", "monorail"));
        config.setNLinkThreshold(16);
        config.setCandidateDistanceMultiplier(1.1);

        return config;
    }

    /**
     * Runs the normal mapping process for the transit schedule and network.
     */
    private static void runNormalMapping() {
        PublicTransitMappingConfigGroup config = createTestPTMConfig();
        TransitSchedule schedule = ScheduleTools.readTransitSchedule(UNMAPPED_SCHEDULE_MERGE);
        Network network = NetworkTools.readNetwork(TEMP_NETWORK_XML.replace(".xml.gz", "_cleaned.xml.gz"));

        PTMapper.mapScheduleToNetwork(schedule, network, config);

        NetworkTools.writeNetwork(network, FINAL_NORMAL_NETWORK_XML);
        ScheduleTools.writeTransitSchedule(schedule, FINAL_NORMAL_SCHEDULE_XML);

        Network2Geojson.run(OUTPUT_COORDINATE_SYSTEM, network, FINAL_NORMAL_NETWORK_GEOJSON);
        Schedule2Geojson.run(OUTPUT_COORDINATE_SYSTEM, schedule, FINAL_NORMAL_SCHEDULE_GEOJSON);

        CheckMappedSchedulePlausibility.run(FINAL_NORMAL_SCHEDULE_XML, FINAL_NORMAL_NETWORK_XML, OUTPUT_COORDINATE_SYSTEM, OUTPUT_DIR + "check/");
    }

    /**
     * Runs the mapping process with shapes for the transit schedule and network.
     */
    private static void runMappingWithShapes() {
//        PublicTransitMappingConfigGroup config = createTestPTMConfig();
//        TransitSchedule schedule = ScheduleTools.readTransitSchedule(UNMAPPED_SCHEDULE_MERGE);
//        Network network = NetworkTools.readNetwork(TEMP_NETWORK_XML.replace(".xml.gz", "_cleaned.xml.gz"));
//
//        ScheduleRoutersFactory routersFactory = new ScheduleRoutersGeojsonShapes.Factory(schedule, network,
//                ShapeTools.readGeojsonShapesFiles(BASE_DIR + "resources/routes.geojson", OUTPUT_COORDINATE_SYSTEM),
//                config.getTransportModeAssignment(), config.getTravelCostType(), 50, 200);
//
//        PTMapper ptMapper = new PTMapper(schedule, network);
//        ptMapper.run(config, null, routersFactory);
//
//        NetworkTools.writeNetwork(network, FINAL_SHAPES_NETWORK_XML);
//        ScheduleTools.writeTransitSchedule(ptMapper.getSchedule(), FINAL_SHAPES_SCHEDULE_XML);
//
//        Network2Geojson.run(OUTPUT_COORDINATE_SYSTEM, network, FINAL_SHAPES_NETWORK_GEOJSON);
//        Schedule2Geojson.run(OUTPUT_COORDINATE_SYSTEM, ptMapper.getSchedule(), FINAL_SHAPES_SCHEDULE_GEOJSON);
//
//        CheckMappedSchedulePlausibility.run(FINAL_SHAPES_SCHEDULE_XML, FINAL_SHAPES_NETWORK_XML, OUTPUT_COORDINATE_SYSTEM, OUTPUT_DIR + "check/");
    }

    

    /**
     * Runs the mapping analysis for the normal mapping.
     */
    public static void mappingAnalysisNormal() {
        runNormalMapping();

//        MappingAnalysis analysis = new MappingAnalysis(
//                ScheduleTools.readTransitSchedule(FINAL_NORMAL_SCHEDULE_XML),
//                NetworkTools.readNetwork(FINAL_NORMAL_NETWORK_XML),
//                ShapeTools.readGeojsonShapesFiles(BASE_DIR + "resources/routes.geojson", OUTPUT_COORDINATE_SYSTEM)
//        );

//        analysis.run();
//        analysis.writeQuantileDistancesCsv(OUTPUT_DIR + "Normal_DistancesQuantile.csv");
    }

    /**
     * Runs the mapping analysis for the mapping with shapes.
     */
    public static void mappingAnalysisWithShapes() {
//        runMappingWithShapes();

//        MappingAnalysis analysis = new MappingAnalysis(
//                ScheduleTools.readTransitSchedule(FINAL_SHAPES_SCHEDULE_XML),
//                NetworkTools.readNetwork(FINAL_SHAPES_NETWORK_XML),
//                ShapeTools.readGeojsonShapesFiles(BASE_DIR + "resources/routes.geojson", OUTPUT_COORDINATE_SYSTEM)
//        );

//        analysis.run();
//        analysis.writeQuantileDistancesCsv(OUTPUT_DIR + "Shapes_DistancesQuantile.csv");
    }

    /**
     * Cuts off the network outside the specified boundary.
     *
     * @param network The network to be cut off.
     * @param boundaryShpFile The boundary shapefile.
     * @throws IOException If there is an I/O error.
     */
    private static void cutOffNetwork(Network network, String boundaryShpFile) throws IOException {
        GeometryFactory geometryFactory = new GeometryFactory();
        STRtree spatialIndex = new STRtree();

        if (boundaryShpFile != null) {
            List<Polygon> polygons = shpPolygonReader.readPolygonsFromFile(boundaryShpFile);
            for (Polygon polygon : polygons) {
                spatialIndex.insert(polygon.getEnvelopeInternal(), polygon);
            }
        }

        List<Link> linksToRemove = new ArrayList<>();
        for (Link link : network.getLinks().values()) {
            if (shouldRemoveLink(link, boundaryShpFile, spatialIndex, geometryFactory)) {
                linksToRemove.add(link);
            }
        }

        for (Link link : linksToRemove) {
            network.removeLink(link.getId());
        }
    }

    /**
     * Determines if a link should be removed based on the boundary shapefile.
     *
     * @param link The link to be checked.
     * @param boundaryShpFile The boundary shapefile.
     * @param spatialIndex The spatial index.
     * @param geometryFactory The geometry factory.
     * @return true if the link should be removed, false otherwise.
     */
    private static boolean shouldRemoveLink(Link link, String boundaryShpFile, STRtree spatialIndex, GeometryFactory geometryFactory) {
        if (link.getAttributes().getAttribute("osm:way:access") != null) {
            String access = (String) link.getAttributes().getAttribute("osm:way:access");
            String bus = (String) link.getAttributes().getAttribute("osm:way:bus");
            if ((access.equals("no") || access.equals("private")) && bus == null) {
                return true;
            }
        }

        if (boundaryShpFile != null) {
            Coord fromCoord = link.getFromNode().getCoord();
            Coord toCoord = link.getToNode().getCoord();
            LineString linkLineString = geometryFactory.createLineString(new Coordinate[]{
                    new Coordinate(fromCoord.getX(), fromCoord.getY()),
                    new Coordinate(toCoord.getX(), toCoord.getY())
            });

            @SuppressWarnings("unchecked")
            List<Polygon> potentialPolygons = spatialIndex.query(linkLineString.getEnvelopeInternal());
            for (Polygon polygon : potentialPolygons) {
                if (polygon.intersects(linkLineString) || polygon.contains(linkLineString)) {
                    return false;
                }
            }
            return true;
        }
        return false;
    }



    /**
     * Reads polygons from a shapefile.
     */
    private static class shpPolygonReader {
        /**
         * Reads polygons from the specified shapefile.
         *
         * @param filePath The path to the shapefile.
         * @return A list of polygons.
         * @throws IOException If there is an I/O error.
         */
        public static List<Polygon> readPolygonsFromFile(String filePath) throws IOException {
            File file = new File(filePath);
            FileDataStore store = FileDataStoreFinder.getDataStore(file);
            SimpleFeatureSource featureSource = store.getFeatureSource();
            FeatureIterator<SimpleFeature> iterator = featureSource.getFeatures().features();
            List<Polygon> polygons = new ArrayList<>();
            try {
                while (iterator.hasNext()) {
                    SimpleFeature feature = iterator.next();
                    Geometry geometry = (Geometry) feature.getDefaultGeometry();
                    if (geometry instanceof MultiPolygon) {
                        for (int i = 0; i < geometry.getNumGeometries(); i++) {
                            polygons.add((Polygon) geometry.getGeometryN(i));
                        }
                    } else if (geometry instanceof Polygon) {
                        polygons.add((Polygon) geometry);
                    }
                }
            } finally {
                iterator.close();
                store.dispose();
            }
            return polygons;
        }
    }
}