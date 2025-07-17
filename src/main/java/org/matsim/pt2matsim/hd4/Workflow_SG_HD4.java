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
import org.matsim.core.network.NetworkUtils;
import org.matsim.pt2matsim.config.OsmConverterConfigGroup;
import org.matsim.pt2matsim.run.Osm2MultimodalNetwork;
import org.matsim.pt2matsim.run.gis.Network2Geojson;
import org.matsim.pt2matsim.tools.NetworkTools;
import org.matsim.utils.objectattributes.attributable.Attributes;
import picocli.CommandLine;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Workflow_SG class handles the workflow for processing and mapping public transit data in Singapore.
 */
@CommandLine.Command(name = "Workflow_SG", description = "Workflow for processing and mapping public transit data in Singapore.")
public class Workflow_SG_HD4 implements MATSimAppCommand {
    private static final Logger log = LogManager.getLogger(Workflow_SG_HD4.class);

    // Directory paths and file names
    @CommandLine.Option(names = "--base-dir", description = "Base directory path", defaultValue = "")
    public static String BASE_DIR;
    
    @CommandLine.Option(names = "--output-dir", description = "Output directory path", defaultValue = "output/")
    private static String OUTPUT_DIR;

    @CommandLine.Option(names = "--temp-network-xml", description = "Temporary network XML file path", defaultValue = "islandwide/walking_cycling_network_epsg3414_singapore_2024_all_accessible_unconnected.xml.gz")
    private static String TEMP_NETWORK_XML;
    
    @CommandLine.Option(names = "--input-osm-file", description = "Input OSM file path", defaultValue ="input/singapore_all_250101_fixed.osm")
    private String INPUT_OSM_FILE;
    
    @CommandLine.Option(names = "--final-normal-network-xml", description = "Final normal network XML file path", defaultValue = "clementi/walking_cycling_network_epsg3414_clementi_2024_public_accessible_unconnected.xml.gz")
    public static String FINAL_NORMAL_NETWORK_XML;

    @CommandLine.Option(names = "--boundary-shp-file", description = "Boundary shapefile path", defaultValue ="input/clementi_PA.shp")
    private String BOUNDARY_SHP_FILE;
    
    @CommandLine.Option(names = "--final-normal-network-geojson", description = "Final normal network GeoJSON file path", defaultValue = "clementi/walking_cycling_network_epsg3414_clementi_2024_public_accessible_unconnected.geojson")
    private static String FINAL_NORMAL_NETWORK_GEOJSON;

    @CommandLine.Option(names = "--output-detailed-link-geometry-csv-file", description = "Output detailed link geometry file path", defaultValue = "islandwide/walking_cycling_network_epsg3414_singapore_2024_all_accessible_unconnected.csv")
    private String OUTPUT_DETAILED_LINK_GEOMETRY_FILE;
    
    @CommandLine.Option(names = "--output-coordinate-system", description = "Output coordinate system", defaultValue = "EPSG:3414")
    private static String OUTPUT_COORDINATE_SYSTEM;


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
        TEMP_NETWORK_XML = OUTPUT_DIR + TEMP_NETWORK_XML;
        INPUT_OSM_FILE = BASE_DIR + INPUT_OSM_FILE;
        FINAL_NORMAL_NETWORK_XML = OUTPUT_DIR + FINAL_NORMAL_NETWORK_XML;
        BOUNDARY_SHP_FILE = BASE_DIR + BOUNDARY_SHP_FILE;
        FINAL_NORMAL_NETWORK_GEOJSON = OUTPUT_DIR + FINAL_NORMAL_NETWORK_GEOJSON;
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

        // Ensure directory exists for temporary network file
        ensureDirectoryExists(TEMP_NETWORK_XML);

        Osm2MultimodalNetwork.run(osmConfig);

        log.info("==============================================================");
        log.info("=      Cutting off network outside boundary of Singapore     =");
        log.info("==============================================================");

        Network network = NetworkTools.readNetwork(osmConfig.getOutputNetworkFile());
        // Adjust network to ensure link speeds are set correctly for walking and cycling
        adjustLinkSpeeds(network, 5.0);


        // Ensure directory exists for geojson file
        ensureDirectoryExists(osmConfig.getOutputNetworkFile().replace(".xml.gz", ".geojson"));
        Network2Geojson.run(OUTPUT_COORDINATE_SYSTEM, network, osmConfig.getOutputNetworkFile().replace(".xml.gz", ".geojson"));

        Network connectedNetwork = NetworkTools.readNetwork(osmConfig.getOutputNetworkFile());
        NetworkUtils.cleanNetwork(connectedNetwork, Set.of("bike", "walk"));
        NetworkUtils.writeNetwork(connectedNetwork, osmConfig.getOutputNetworkFile().replace("unconnected", "connected"));
        Network2Geojson.run(OUTPUT_COORDINATE_SYSTEM, connectedNetwork, osmConfig.getOutputNetworkFile().replace("unconnected", "connected").replace(".xml.gz", ".geojson"));

        cutOffNetwork(network, BOUNDARY_SHP_FILE);

        // Ensure directory exists for final network file
        ensureDirectoryExists(FINAL_NORMAL_NETWORK_XML);
        NetworkUtils.writeNetwork(network, FINAL_NORMAL_NETWORK_XML);

        log.info("==> Network cut off successfully!");
        log.info("==============================================================");

        // Ensure directory exists for final geojson file
        ensureDirectoryExists(FINAL_NORMAL_NETWORK_GEOJSON);
        Network2Geojson.run(OUTPUT_COORDINATE_SYSTEM, network, FINAL_NORMAL_NETWORK_GEOJSON);

        //Clean network to make it connected
        log.info("==============================================================");
        log.info("=      Cleaning network to make it connected                  =");
        log.info("==============================================================");
        NetworkUtils.cleanNetwork(network, Set.of("bike", "walk"));
        NetworkUtils.writeNetwork(network, FINAL_NORMAL_NETWORK_XML.replace("unconnected", "connected"));
        Network2Geojson.run(OUTPUT_COORDINATE_SYSTEM, network, FINAL_NORMAL_NETWORK_GEOJSON.replace("unconnected", "connected"));
        log.info("==> Network cleaned successfully!");
        log.info("===========================================================");

        return 0;
    }

    private void adjustLinkSpeeds(Network network, double maxWalkingSpeed_kmh) {
        AtomicInteger adjustedCount = new AtomicInteger();
        for (Link link : network.getLinks().values()) {
            Attributes attrs = link.getAttributes();
            attrs.getAsMap().forEach((key, value) -> {
                // Check if it's a footway on a bridge (with null check)
                if (key.equals("osm:way:highway") && value.equals("footway")) {
                    Object bridgeAttr = attrs.getAttribute("osm:way:bridge");
                    if (bridgeAttr != null && bridgeAttr.equals("yes")) {
                        // Adjust cycling speed for cycleways on bridges
                        link.setFreespeed(maxWalkingSpeed_kmh / 3.6); // Convert km/h to m/s
//                        log.info("Adjusted freespeed for cycleway bridge link: {}", link.getId());
                        adjustedCount.getAndIncrement();
                    }
                }

                // Handle maxspeed attributes
                if (key.contains("maxspeed") && value instanceof Number) {
                    if (key.contains("official")) {
                        // Adjust cycling speed for links with maxspeed
                        double maxSpeed = ((Number) value).doubleValue();
                        link.setFreespeed(maxSpeed / 3.6); // Convert km/h to m/s
                        adjustedCount.getAndIncrement();
                    } else if (key.contains("advisory")) {
                        // Adjust walking speed for links with maxspeed
                        double maxSpeed = ((Number) value).doubleValue();
                        link.setFreespeed(maxSpeed / 3.6); // Convert km/h to m/s
                    }
                }
            });
        }
        log.info("Adjusted link speeds for walking and cycling on bridges: {}", adjustedCount.get());
    }

    /**
     * Ensures that the directory for a file exists, creating it if necessary
     *
     * @param filePath Path to the file
     */
    private void ensureDirectoryExists(String filePath) {
        try {
            Path path = Paths.get(filePath);
            if (path.getParent() != null) {
                Files.createDirectories(path.getParent());
                log.info("Created directory structure for: {}", path.getParent());
            }
        } catch (IOException e) {
            log.error("Failed to create directory for: {}", filePath, e);
        }
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