/* Acts as constructor for Phenotype objects */
/* A completely static class */

package org.antigen.phenotype;

import org.antigen.core.Parameters;

public class PhenotypeFactory {

  public static String GEOMETRIC = "geometric";
  public static String GEOMETRIC3D = "geometric3d";
  public static String GEOMETRIC10D = "geometric10d";
  public static String SEQUENCE = "sequence";
  public static String GEOMETRIC_SEQ = "geometricSeq";

  // returns newly instantiated Phenotype objects of type according to Parameters.phenotypeSpace
  public static Phenotype makeVirusPhenotype() {

    Phenotype p = null;
    if (GEOMETRIC.equals(Parameters.phenotypeSpace)) {
      p = new GeometricPhenotype();
    }
    if (GEOMETRIC3D.equals(Parameters.phenotypeSpace)) {
      p = new GeometricPhenotype3D();
    }
    if (GEOMETRIC10D.equals(Parameters.phenotypeSpace)) {
      p = new GeometricPhenotype10D();
    }
    if (GEOMETRIC_SEQ.equals(Parameters.phenotypeSpace)) {
      p = new GeometricSeqPhenotype();
    }
    return p;
  }

  // returns newly instantiated Phenotype objects of type according to Parameters.phenotypeSpace
  public static Phenotype makeHostPhenotype() {

    Phenotype p = null;
    if (GEOMETRIC.equals(Parameters.phenotypeSpace)) {
      p = new GeometricPhenotype(Parameters.initialTraitA, 0);
    }
    if (GEOMETRIC3D.equals(Parameters.phenotypeSpace)) {
      p = new GeometricPhenotype3D(Parameters.initialTraitA, 0, 0);
    }
    if (GEOMETRIC10D.equals(Parameters.phenotypeSpace)) {
      double[] traits = {Parameters.initialTraitA, 0, 0, 0, 0, 0, 0, 0, 0, 0};
      p = new GeometricPhenotype10D(traits);
    }
    if (GEOMETRIC_SEQ.equals(Parameters.phenotypeSpace)) {
      String startingSequence = Parameters.startingSequence;
      if (startingSequence == null) {
        p = new GeometricSeqPhenotype(Parameters.initialTraitA, 0);
      } else {
        p =
            new GeometricSeqPhenotype(
                Parameters.initialTraitA, 0, Parameters.startingSequence.toCharArray());
      }
    }
    return p;
  }

  // returns newly instantiated Phenotype objects of type according to Parameters.phenotypeSpace
  public static Phenotype makeArbitaryPhenotype(double x, double y) {

    Phenotype p = null;
    if (GEOMETRIC.equals(Parameters.phenotypeSpace)) {
      p = new GeometricPhenotype(x, y);
    }
    if (GEOMETRIC_SEQ.equals(Parameters.phenotypeSpace)) {
      p = new GeometricSeqPhenotype(x, y);
    }
    return p;
  }
}
