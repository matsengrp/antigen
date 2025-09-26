package org.antigen.host;

import org.antigen.core.Parameters;
import org.antigen.phenotype.GeometricPhenotype;
import org.antigen.phenotype.GeometricSeqPhenotype;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

public class TestHostMostRecentInfection {
    
    @Before
    public void setUp() {
        // Disable initial immunity to ensure hosts start naive
        Parameters.initialPrR = 0.0;
    }
    
    @Test
    public void testNaiveHostCentroid() {
        Host host = new Host();
        double[] centroid = host.getImmunityCoordinatesCentroid();
        assertNull("Naive host should return null centroid", centroid);
    }
    
    @Test
    public void testSingleInfectionMostRecent() {
        Host host = new Host();
        GeometricPhenotype phenotype = new GeometricPhenotype(2.0, 3.0);
        host.addToHistory(phenotype);
        
        double[] mostRecent = host.getImmunityCoordinatesCentroid();
        assertNotNull(mostRecent);
        assertEquals(2.0, mostRecent[0], 1e-10);
        assertEquals(3.0, mostRecent[1], 1e-10);
    }
    
    @Test
    public void testMultipleInfectionsMostRecent() {
        Host host = new Host();
        
        // Add first infection
        GeometricPhenotype p1 = new GeometricPhenotype(1.0, 2.0);
        host.addToHistory(p1);
        
        // Add second infection (more recent)
        GeometricPhenotype p2 = new GeometricPhenotype(3.0, 4.0);
        host.addToHistory(p2);
        
        double[] mostRecent = host.getImmunityCoordinatesCentroid();
        assertNotNull(mostRecent);
        assertEquals(3.0, mostRecent[0], 1e-10); // Should return coords of p2 (most recent)
        assertEquals(4.0, mostRecent[1], 1e-10); // Should return coords of p2 (most recent)
    }
    
    @Test
    public void testGeometricSeqPhenotypeCentroid() {
        Host host = new Host();
        GeometricSeqPhenotype phenotype = new GeometricSeqPhenotype(5.0, 7.0);
        host.addToHistory(phenotype);
        
        double[] centroid = host.getImmunityCoordinatesCentroid();
        assertNotNull(centroid);
        assertEquals(5.0, centroid[0], 1e-10);
        assertEquals(7.0, centroid[1], 1e-10);
    }
    
    @Test
    public void testMixedPhenotypesMostRecent() {
        Host host = new Host();
        
        // Add GeometricPhenotype first
        GeometricPhenotype gp = new GeometricPhenotype(2.0, 4.0);
        host.addToHistory(gp);
        
        // Add GeometricSeqPhenotype (more recent)
        GeometricSeqPhenotype gsp = new GeometricSeqPhenotype(6.0, 8.0);
        host.addToHistory(gsp);
        
        double[] mostRecent = host.getImmunityCoordinatesCentroid();
        assertNotNull(mostRecent);
        assertEquals(6.0, mostRecent[0], 1e-10); // Should return coords of gsp (most recent)
        assertEquals(8.0, mostRecent[1], 1e-10); // Should return coords of gsp (most recent)
    }
    
    @Test
    public void testMostRecentOrderMatters() {
        Host host = new Host();
        
        // Add three infections in sequence
        GeometricPhenotype p1 = new GeometricPhenotype(1.0, 1.0);
        GeometricPhenotype p2 = new GeometricPhenotype(5.0, 5.0);
        GeometricPhenotype p3 = new GeometricPhenotype(9.0, 9.0);
        
        host.addToHistory(p1);
        host.addToHistory(p2);
        host.addToHistory(p3);
        
        double[] mostRecent = host.getImmunityCoordinatesCentroid();
        assertNotNull(mostRecent);
        assertEquals(9.0, mostRecent[0], 1e-10); // Should return coords of p3 (last added)
        assertEquals(9.0, mostRecent[1], 1e-10); // Should return coords of p3 (last added)
    }
}