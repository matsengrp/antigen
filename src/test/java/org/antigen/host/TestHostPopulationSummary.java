package org.antigen.host;

import org.antigen.core.Parameters;
import org.antigen.phenotype.GeometricPhenotype;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

public class TestHostPopulationSummary {
    
    private HostPopulation population;
    
    @Before
    public void setUp() {
        Parameters.demeCount = 1;
        Parameters.initialNs = new int[]{100};
        Parameters.initialPrR = 0.0; // Disable initial immunity
        population = new HostPopulation(0);
    }
    
    @Test
    public void testAllNaivePopulation() {
        ImmunitySummary summary = population.getPopulationImmunitySummary(10);
        
        assertEquals("All hosts should be naive", 1.0, summary.getNaiveFraction(), 1e-10);
        assertEquals("No experienced hosts", 0, summary.getExperiencedHosts());
        assertEquals("Total sampled should be 10", 10, summary.getTotalSampled());
        assertFalse("Should not have valid centroid", summary.hasValidCentroid());
        assertTrue("Centroid should be NaN", Double.isNaN(summary.getCentroid()[0]));
        assertTrue("Centroid should be NaN", Double.isNaN(summary.getCentroid()[1]));
    }
    
    @Test
    public void testMixedPopulation() {
        // Add some infected hosts with known immunity centroids
        for (int i = 0; i < 5; i++) {
            Host host = population.getRandomHost();
            GeometricPhenotype p = new GeometricPhenotype(i * 1.0, i * 2.0);
            host.addToHistory(p);
        }
        
        // Sample the population
        ImmunitySummary summary = population.getPopulationImmunitySummary(20);
        
        // Should have some naive and some experienced hosts
        assertTrue("Should have some naive hosts", summary.getNaiveFraction() > 0.0);
        assertTrue("Should have some experienced hosts", summary.getNaiveFraction() < 1.0);
        assertTrue("Should have experienced hosts", summary.getExperiencedHosts() > 0);
        assertEquals("Total sampled should be 20", 20, summary.getTotalSampled());
        
        if (summary.hasValidCentroid()) {
            assertFalse("Centroid should not be NaN", Double.isNaN(summary.getCentroid()[0]));
            assertFalse("Centroid should not be NaN", Double.isNaN(summary.getCentroid()[1]));
        }
    }
    
    @Test
    public void testImmunitySummaryProperties() {
        ImmunitySummary summary = new ImmunitySummary(
            new double[]{1.0, 2.0}, 
            0.3, 
            100, 
            70
        );
        
        assertEquals("Centroid A coordinate", 1.0, summary.getCentroid()[0], 1e-10);
        assertEquals("Centroid B coordinate", 2.0, summary.getCentroid()[1], 1e-10);
        assertEquals("Naive fraction", 0.3, summary.getNaiveFraction(), 1e-10);
        assertEquals("Total sampled", 100, summary.getTotalSampled());
        assertEquals("Experienced hosts", 70, summary.getExperiencedHosts());
        assertTrue("Should have valid centroid", summary.hasValidCentroid());
    }
    
    @Test
    public void testNaNHandling() {
        ImmunitySummary summary = new ImmunitySummary(
            new double[]{Double.NaN, Double.NaN}, 
            1.0, 
            50, 
            0
        );
        
        assertFalse("Should not have valid centroid", summary.hasValidCentroid());
        assertTrue("Centroid should be NaN", Double.isNaN(summary.getCentroid()[0]));
        assertTrue("Centroid should be NaN", Double.isNaN(summary.getCentroid()[1]));
    }
}