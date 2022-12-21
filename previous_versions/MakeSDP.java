/**
 * This program creates the SDP used for proving Property (T) for Aut(F_4).
 * @author Martin Nitsche
 * @version LapSym-MaxEq-01
 * The mathematical details are explained in the accompanying article and the sources cited therein.
 * 
 * This program will run for a few seconds and create multiple files of combined size ca. 200mb:
 *   - "A.dat", "b.dat", "c.dat" describe the SDP problem in the usual (sparse matrix) format.
 *     Technically, this data would suffice for solving the SDP problem and proving Property (T).
 *   - "diagonalAction.dat", "permutationAction.dat", "blockSizes.dat" contain information
 *     about the the symmetry of the problem, used by the MATLAB script to speed up the calculation.
 *   - "points.txt", "distances.txt", "blocks.txt" are meant for debugging or experimentation.
 *     They document which group elements correspond to which matrix/vector entries in A, b, c.
 * 
 * To run this code, install the Java SDK, compile with "javac MakeSDP.java",
 * create a folder "data" in the same directory, and run with "java MakeSDP".
 * To solve the SDP use the accompanying MATLAB script.
 * 
 * We tried to keep the code as self-contained as possible and to use little abstraction or Java-
 * specific features. Documentation for the Java library can be found at docs.oracle.com, e.g. 
 * https://docs.oracle.com/en/java/javase/13/docs/api/java.base/java/util/Hashtable.html
 * 
 * If you want to experiment with this code, you can start by changing the blocksToInclude variable
 * to contain a slightly different subset of the Strings in "blocks.txt".
 * But be aware that even small changes can make the SDP significantly more complex.
 * In order to work with non-symmetric sets of points, larger adjustments will have to be made.
 */

import java.util.Arrays;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.ArrayList;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

// The main class, producing the SDP.
// Note: Although we split the code up into functions for readability, MakeSDP is really just a linear script.
class MakeSDP {

  // The directory where the output is written
  String outputDirectory = "./data/";

  // The set of points on which the semidefinite variable in the SDP problem lives consists of equivalence classes
  // of automorphisms. We identify these classes with fixed representatives obtained from pointClassRepresentative.
  // The order of the points will remain fixed and will later determine the order of the entries of A and c.
  // The identity automorphism will not be included, but as the origin it will always be part of the SDP problem.
  // points will contain the list of points, pointNumbers will contain for each point the position in the list.
  Automorphism[]                   points;       // array for fixed order and fast random access
  Hashtable<Automorphism, Integer> pointNumbers; // Hashtable for fast look-up

  // The list of points will be a subset of a list of point candidates. A point equivanece class shall be a point
  // candidate if it intersects the 3-ball around identity.
  // For each point equivalence class in pointCandidates we also store the inverse of some automorphism in the
  // equivalence class. These inverses are needed to compute the distances. It is faster and easier to store them
  // than to compute them on the fly later.
  int POINT_CANDIDATES_RADIUS = 3;
  ArrayList<Automorphism> pointCandidates;             // ArrayList for fixed order and variable size
  Hashtable<Automorphism, Automorphism> inverseLookup; // Hashtable for fast look-up

  // The distances between the points, from which the SDP objective and constraints are calculated, are also given
  // by equivalence classes of automorphisms.  If two points are represented by automorphisms a, b, then their distance
  // is determined by the equivalence class of a^{-1}b, which does not depend on the choice of the representatives.
  // We identify the distance equivalence classes with fixed representatives obtained from distanceClassRepresentative.
  // We will collect them as we encounter them and number them. The classes for identity and the generators are special
  // and will be given negative numbers. The other numbers determine the order of the rows in A and b.
  // All distance classes, except those of identity and the generators, are also stored in order in the distances list.
  int ZERO_DISTANCE_NUMBER = -2;
  int UNIT_DISTANCE_NUMBER = -1;
  Hashtable<Automorphism, Integer> distanceNumbers; // Hashtable for fast look-up
  ArrayList<Automorphism>          distances;       // ArrayList for fixed order and variable size

  // The distance equivalence relation is coarser than the point equivalence relation and invariant under the symmetries.
  // We use it to partition the points obtained from the 3-ball around the identity automorphism into blocks.
  // The blocks represented by the listed automorphisms will be included in the set of points, blockwise.
  // For human readability, we do not work with the Automorphism instances but their unique String representation.
  Hashtable<String, ArrayList<Automorphism>> pointBlocks;
  String[] blocksToInclude = {
    "ab, a, c, d, "  , "ab, cd, a, c, " , "ab, ac, a, d, " , "ab, Ac, a, d, " , "abA, a, c, d, " , // from 2-ball
    "abc, a, b, d, " , "aba, a, c, d, " , "aab, a, c, d, " , "abc, a, c, d, " , "abc, ab, a, d, ",
    "abc, ab, b, d, ", "abA, ab, c, d, ", "aab, ab, c, d, ",
    "ab, ac, ad, a, ", "ab, ac, Ad, a, ", "abA, ac, a, d, ", "aab, ac, a, d, "     // selected blocks from 3-ball
  };

  // To speed up the computation we make use of the diagonal and the permutation subgroup of Aut(F_4), which
  // together generate the symmetry group of the SDP problem. To save time we precompute the action.
  // The i-th entry in Automorphism.diagonalElements sends the j-th entry in points to the point at position
  // diagonalAction[i][j], and similarly for the permutation subgroup.
  // We also precompute the distance numbers for the distance from identity to each point and store them.
  int[][] diagonalAction;
  int[][] permutationAction;
  int[]   distanceToPoint;


  // The main method, called when the program is executed.
  // It simply executes once each function below, in the order that they are written.
  public static void main(String[] args) {
    MakeSDP sdpMaker = new MakeSDP();
    System.out.print("Creating 3-ball... ");
    sdpMaker.createPointCandidates();
    System.out.print("" + sdpMaker.pointCandidates.size() + " point candidates\n");
    System.out.print("Selecting points... ");
    sdpMaker.filterPointsByBlocks();
    System.out.print("" + sdpMaker.points.length + " points selected\n");
    System.out.print("Making precomputations (symmetry of points, distances from origin)... ");
    sdpMaker.makePrecomputations();
    System.out.print("done\n");
    try {
      System.out.print("Computing objective and constraints... ");
      sdpMaker.computeAndWriteAc();
      System.out.print("" + sdpMaker.distances.size() + " constraints\n");
      System.out.print("Writing constraints right hand side\n");
      sdpMaker.computeAndWriteb();
      System.out.print("Writing symmetry data\n");
      sdpMaker.writeSymmetryData();
      System.out.print("Writing documentation data\n");
      sdpMaker.writeDocumentation(); }
    catch(IOException e) { System.err.println("caught: " + e); } // catch possible errors from writing to disk
  }


  // Called by "new MakeSDP()", initializes two hashtables.
  // (The other arrays, lists and hashtables will be initialized by the functions dedicated to filling them.)
  public MakeSDP() {
    inverseLookup = new Hashtable<Automorphism, Automorphism>();
    distanceNumbers = new Hashtable<Automorphism, Integer>();
    Automorphism identity = Automorphism.identity();
    Automorphism someGenerator = Automorphism.elementaryAutomorphism(1, 2);
    distanceNumbers.put(Automorphism.distanceClassRepresentative(identity     , identity     ), ZERO_DISTANCE_NUMBER);
    distanceNumbers.put(Automorphism.distanceClassRepresentative(someGenerator, someGenerator), UNIT_DISTANCE_NUMBER);
  }


  // Fills pointCandidates with all those point equivalence classes that intersect the ball of radius
  // POINT_CANDIDATES_RADIUS around identity (with identity itself excluded).
  // The order of the elements in the resulting array will always be the same.
  // For each point candidate we also store an inverse in inverseLookup.
  void createPointCandidates() {

    // We begin by making a list of the generators for SAut(F_4) and their inverses.
    // We do this by hand, the Automorphism class does not have a method for finding inverses.
    ArrayList<Automorphism> generators = new ArrayList<Automorphism>(); // the elementary automorphisms generating SAut(F_4)
    ArrayList<Automorphism> generatorInverses = new ArrayList<Automorphism>(); // inverses of generators in correct order
    for (int i = 1; i <= 4; i++) {
      for (int j = 1; j <= 4; j++) {
        if (i == j) continue; // no elementary automorphisms for i=j, skip iteration
        generators.add(Automorphism.elementaryAutomorphism( i,  j));
        generatorInverses.add(Automorphism.elementaryAutomorphism( i, -j));
        generators.add(Automorphism.elementaryAutomorphism( i, -j));
        generatorInverses.add(Automorphism.elementaryAutomorphism( i,  j));
        generators.add(Automorphism.elementaryAutomorphism(-i,  j));
        generatorInverses.add(Automorphism.elementaryAutomorphism(-i, -j));
        generators.add(Automorphism.elementaryAutomorphism(-i, -j));
        generatorInverses.add(Automorphism.elementaryAutomorphism(-i,  j));
      }
    }

    // Now we create the candidates, first from the 1-ball, then increasing the radius step by step.
    pointCandidates = new ArrayList<Automorphism>();
    HashSet<Automorphism> candidateSet = new HashSet<Automorphism>();     // fast look-up, which candidates we already met
    for (int i = 0; i < generators.size(); i++) {
      Automorphism pointRepresentative = generators.get(i).pointClassRepresentative();
      if (candidateSet.contains(pointRepresentative)) continue;           // skip iteration if we already have this point
      pointCandidates.add(pointRepresentative);
      candidateSet.add(pointRepresentative);
      inverseLookup.put(pointRepresentative, generatorInverses.get(i));
    }
    for (int i = 2; i <= POINT_CANDIDATES_RADIUS; i++) {
      ArrayList<Automorphism> nextLayer = new ArrayList<Automorphism>();  // new points to add
      for (Automorphism p : pointCandidates) {                            // for every point of the (i-1)-ball (in order)...
        for (int j = 0; j < generators.size(); j++) {                     // ...we multiply with every generator
          // the multiplication order guarantees that the new point does not depend on the representative for point p
          Automorphism prod = Automorphism.product(generators.get(j), p); // multiplication order is important, see below
          Automorphism pointRepresentative = prod.pointClassRepresentative();
          if (candidateSet.contains(pointRepresentative)) continue;       // skip iteration if we already have this point
          nextLayer.add(pointRepresentative);                             // do not add to pointCandidates while inside middle loop
          candidateSet.add(pointRepresentative);
          Automorphism anInverse = Automorphism.product(inverseLookup.get(p), generatorInverses.get(j));
          inverseLookup.put(pointRepresentative, anInverse);              // correct because of multiplication order, see above
        }
      }
      pointCandidates.addAll(nextLayer);
    }
  }


  // Fills pointBlocks with pointCandidates sorted into blocks, according to the distance equivalence relation.
  // Also fills points with the pointCandidates from blocksToInclude (without changing the order).
  void filterPointsByBlocks() {
    pointBlocks = new Hashtable<String, ArrayList<Automorphism>>();
    for (Automorphism candidate : pointCandidates) { // go through pointCandidates, in order
      String blockString = Automorphism.distanceClassRepresentative(candidate, inverseLookup.get(candidate)).toString();
      ArrayList<Automorphism> blockElements = pointBlocks.get(blockString);
      if (blockElements == null) {                   // if blockSizes does not yet have an entry for s, create one
        blockElements = new ArrayList<Automorphism>();
        pointBlocks.put(blockString, blockElements); }
      blockElements.add(candidate);
    }
    int numberOfPoints = 0;                          // total number of points to include
    for (int i = 0; i < blocksToInclude.length; i++) {
      numberOfPoints += pointBlocks.get(blocksToInclude[i]).size(); }
    points = new Automorphism[numberOfPoints];
    pointNumbers = new Hashtable<Automorphism, Integer>();
    int nextNumber = 0;                              // number of next point to be added
    for (int i = 0; i < blocksToInclude.length; i++) {
      ArrayList<Automorphism> blockElements = pointBlocks.get(blocksToInclude[i]);
      for (int j = 0; j < blockElements.size(); j++) {
        points[nextNumber] = blockElements.get(j);
        pointNumbers.put(blockElements.get(j), nextNumber);
        nextNumber++;
      }
    }
  }


  // Fills diagonalAction and permutationAction with the action of the diagonal, respectively permutation subgroups
  // on the points.
  // Also computes the distances from all points to identity and stores the distance numbers in order in
  // distanceToPoint. When new distances are encountered, distances and distanceNumbers are updated.
  // (This naive implementation does not use the symmetry to speed up the computation.)
  void makePrecomputations() {
    diagonalAction    = new int[Automorphism.DIAGONALS.length   ][points.length];
    permutationAction = new int[Automorphism.PERMUTATIONS.length][points.length];
    for (int i = 0; i < points.length; i++) {
      Automorphism[] diagOrbit = points[i].diagonalActionOrbit();
      for (int j = 0; j < diagOrbit.length; j++) {
        diagonalAction[j][i] = pointNumbers.get(diagOrbit[j].pointClassRepresentative()); }
      Automorphism[] permOrbit = points[i].permutationActionOrbit();
      for (int j = 0; j < permOrbit.length; j++) {
        permutationAction[j][i] = pointNumbers.get(permOrbit[j].pointClassRepresentative()); }
    }
    distances = new ArrayList<Automorphism>();
    distanceToPoint = new int[points.length];
    for (int i = 0; i < points.length; i++) {
      Automorphism distanceRepresentative = Automorphism.distanceClassRepresentative(points[i], inverseLookup.get(points[i]));
      if (!distanceNumbers.containsKey(distanceRepresentative)) { // if not yet encountered, add distance
        int nextNumber = distances.size();
        distanceNumbers.put(distanceRepresentative, nextNumber);
        distances.add(distanceRepresentative); }
      distanceToPoint[i] = distanceNumbers.get(distanceRepresentative);
    }
  }


  // The key task: Computes the constraint matrix A and the vector c, which comes from the unit distance.
  // We write A and c to disk as we go without keeping them in memory. Thus, we could process even the whole 3-ball
  // without running out of RAM (requires >100Gb disk space).
  // We make use of the symmetry to reduce computing time, with the two arrays pointLocked and partnerLocked keeping
  // track of which entries were dealt with earlier. Note that partnerLocked is a local variable inside the outer
  // loop. This way we avoid having a two-dimensional array of size (points.length*points.length).
  void computeAndWriteAc() throws IOException {
    FileWriter writer_A = new FileWriter(new File(outputDirectory + "A.dat"), false);
    FileWriter writer_c = new FileWriter(new File(outputDirectory + "c.dat"), false);
    boolean[] pointLocked = new boolean[points.length];     // for which points we have already computed distances to all other points
    for (int i = 0; i < points.length; i++) {
      if (pointLocked[i]) continue;                         // if point is locked, skip this iteration
      boolean[] partnerLocked = new boolean[points.length]; // for which points we have already computed the distance to points[i]
      for (int j = 0; j < points.length; j++) {
        if (pointLocked[j] || partnerLocked[j]) continue;

        // The (i,j)-matrix entry of the semidefinite variable should be
        //   sum of distance to points from identity, minus distance between them
        // If the same distance variable appears multiple times in this calculation, we must write it only once in the
        // constraint matrix with the appropriate coefficient.
        // If the distance number ZERO_DISTANCE_NUMBER appears we disregard it.
        // If the distance number UNIT_DISTANCE_NUMBER appears, we will set it apart for the constraint vector.
        Automorphism prod = Automorphism.product(inverseLookup.get(points[i]), points[j]);
        Automorphism inverse = Automorphism.product(inverseLookup.get(points[j]), points[i]);
        Automorphism distanceRepresentative = Automorphism.distanceClassRepresentative(prod, inverse);
        if (!distanceNumbers.containsKey(distanceRepresentative)) { // if not yet encountered, add distance
          int nextNumber = distances.size();
          if (nextNumber % 1000 == 0) {              // give regular feedback to user, proxy for total progress
            System.out.print(nextNumber + " "); }
          distanceNumbers.put(distanceRepresentative, nextNumber);
          distances.add(distanceRepresentative); }
        int numA, numB, numC;                        // the distance numbers needed for the (i,j)-entry
        numA = distanceToPoint[i];
        numB = distanceToPoint[j];
        numC = distanceNumbers.get(distanceRepresentative);
        int valA = 1, valB = 1, valC = -1, valU = 0; // coefficients for the distances, valU is for unit distance
        if (numA == UNIT_DISTANCE_NUMBER) valU += valA;
        if (numB == UNIT_DISTANCE_NUMBER) valU += valB;
        if (numC == UNIT_DISTANCE_NUMBER) valU += valC;
        if (numA == UNIT_DISTANCE_NUMBER || numA == ZERO_DISTANCE_NUMBER) valA = 0;
        if (numB == UNIT_DISTANCE_NUMBER || numB == ZERO_DISTANCE_NUMBER) valB = 0;
        if (numC == UNIT_DISTANCE_NUMBER || numC == ZERO_DISTANCE_NUMBER) valC = 0;
        if (numB == numA) { valA += valB; valB = 0; }
        if (numC == numA) { valA += valC; valC = 0; }
        if (numC == numB) { valB += valC; valC = 0; }

        // Using the symmetry we can write the contraints not only for the entry (i,j) but for its whole G-orbit,
        // and also for the mirror entries, i.e. the G-orbit of (j,i). Because the G-action might not be free,
        // we must keep track of which entries we have already written. Note, however that we do not have to check
        // for pointLocked or partnerLocked, because if one of these flags were set for any point in the G-orbit of
        // points[i] or points[j], they would have been set for points[i], points[j] respectively, which is checked above.
        HashSet<Integer> lockedEntries = new HashSet<Integer>();     // entries already written while iterating over symmetry group G
        for (int k = 0; k < Automorphism.DIAGONALS.length; k++) {
          for (int l = 0; l < Automorphism.PERMUTATIONS.length; l++) {
            int curRow = diagonalAction[k][permutationAction[l][i]]; // image of i-th point under action of current group element
            int curCol = diagonalAction[k][permutationAction[l][j]]; // image of j-th point under action of current group element
            int entry = curRow * points.length + curCol + 1;         // resulting entry in vectorized matrix, +1 for MATLAB indexing
            int mirrorEntry = curCol * points.length + curRow + 1;   // transposed entry in vectorized matrix
            if (lockedEntries.contains(entry)) continue;             // if entry locked skip this iteration
            if (curRow == i) partnerLocked[curCol] = true;           // lock partners, but careful: partnerLocked is only for points[i]
            if (curCol == i) partnerLocked[curRow] = true;
            lockedEntries.add(entry);
            lockedEntries.add(mirrorEntry);
            if (valA != 0) writer_A.write((numA+1) + " " + entry + " " + valA + "\n"); // for a sparse matrix we skip trivial values
            if (valB != 0) writer_A.write((numB+1) + " " + entry + " " + valB + "\n");
            if (valC != 0) writer_A.write((numC+1) + " " + entry + " " + valC + "\n");
            if (valU != 0) writer_c.write(entry + " 1 " + valU + "\n");
            if (curRow == curCol) continue;                          // do not write the diagonal entries twice, skip iteration
            if (valA != 0) writer_A.write((numA+1) + " " + mirrorEntry + " " + valA + "\n");
            if (valB != 0) writer_A.write((numB+1) + " " + mirrorEntry + " " + valB + "\n");
            if (valC != 0) writer_A.write((numC+1) + " " + mirrorEntry + " " + valC + "\n");
            if (valU != 0) writer_c.write(mirrorEntry + " 1 " + valU + "\n");
          }
        }
      }

      // After finishing the inner loop update pointLocked: lock the processed point and its whole G-orbit.
      for (int k = 0; k < Automorphism.DIAGONALS.length; k++) {
        for (int l = 0; l < Automorphism.PERMUTATIONS.length; l++) {
          pointLocked[diagonalAction[k][permutationAction[l][i]]] = true; }
      }
    }

    // if necessary, write a final zero to indicate the size of c, for A this is also always necessary
    if (distanceToPoint[points.length-1] > 0) {
      writer_c.write("" + (points.length*points.length) + " 1 0\n"); }
    writer_A.write("" + distances.size() + " " + (points.length*points.length) + " 0\n");
    writer_A.close();
    writer_c.close();
  }


  // Calculates b, the right hand side of the constraints, and writes it to disk.
  // Because of the symmetry, it suffices to multiply all generators with one fixed generator.
  // We first collect all entries, then output as sparse matrix.
  void computeAndWriteb() throws IOException {
    int[] bEntries = new int[distances.size()]; // entries of the b vector
    Automorphism[] generators = Automorphism.generators();
    for (Automorphism a : generators) {
      Automorphism prod = Automorphism.product(generators[0], a).pointClassRepresentative();
      int distNum = distanceNumbers.get(Automorphism.distanceClassRepresentative(prod, inverseLookup.get(prod)));
      if (distNum != ZERO_DISTANCE_NUMBER && distNum != UNIT_DISTANCE_NUMBER) {
        bEntries[distNum]++; }
    }
    FileWriter writer_b = new FileWriter(new File(outputDirectory + "b.dat"), false);
    for (int i = 0; i < bEntries.length; i++) {
      if (bEntries[i] != 0 || i == bEntries.length-1) {            // in sparse matrix format write only the last zero entry
        writer_b.write("1 " + (i+1) + " " + bEntries[i] + "\n"); } // i+1 because MATLAB starts indexing at 1
    }
    writer_b.close();
  }


  // Outputs information about the symmetry of the SDP problem: The sizes of the invariant blocks in the
  // order that they appear in points, and the action of the permutation and the diagonal subgroup on points.
  // The latter is encoded as one matrix each, where the i-th row is the vectorized permutation matrix for
  // the i-th element.
  void writeSymmetryData() throws IOException {
    FileWriter writerSizes        = new FileWriter(new File(outputDirectory + "blockSizes.dat"       ), false);
    FileWriter writerDiagonals    = new FileWriter(new File(outputDirectory + "diagonalAction.dat"   ), false);
    FileWriter writerPermutations = new FileWriter(new File(outputDirectory + "permutationAction.dat"), false);
    for (int i = 0; i < blocksToInclude.length; i++) {
      writerSizes.write("" + (i+1) + " 1 " + pointBlocks.get(blocksToInclude[i]).size() + "\n"); }
    for (int i = 0; i < diagonalAction.length; i++) {
      for (int j = 0; j < points.length; j++) {
        int entry = j*points.length+diagonalAction[i][j]+1;
        writerDiagonals.write((i+1) + " " + entry + " 1\n"); }
    }
    for (int i = 0; i < permutationAction.length; i++) {
      for (int j = 0; j < points.length; j++) {
        int entry = j*points.length+permutationAction[i][j]+1;
        writerPermutations.write((i+1) + " " + entry + " 1\n"); }
    }
    writerSizes.close();
    writerDiagonals.close();
    writerPermutations.close();
  }


  // Outputs information about the block decomposition of the 3-ball and which blocks we included as points.
  // Also output a list of the point and distance equivalence classes in our fixed order
  // (which determines the order of the rows and columns of A, b and c).
  void writeDocumentation() throws IOException {
    FileWriter writerBlocks    = new FileWriter(new File(outputDirectory + "blocks.txt"   ), false);
    FileWriter writerPoints    = new FileWriter(new File(outputDirectory + "points.txt"   ), false);
    FileWriter writerDistances = new FileWriter(new File(outputDirectory + "distances.txt"), false);
    for (String s : pointBlocks.keySet()) {
      writerBlocks.write(s + "(" + pointBlocks.get(s).size() + " points)");
      if (pointNumbers.containsKey(pointBlocks.get(s).get(0))) {
        writerBlocks.write(" included\n"); }
      else {
        writerBlocks.write("\n"); }
    }
    for (int i = 0; i < points.length; i++) {
      writerPoints.write("" + (i+1) + " # " + points[i] + "\n"); }
    writerPoints.close();
    for (int i = 0; i < distances.size(); i++) {
      writerDistances.write("" + (i+1) + " # " + distances.get(i) + "\n"); }
    writerBlocks.close();
    writerPoints.close();
    writerDistances.close();
  }

}


// The elements of Aut(F_4) are represented by instances of the Automorphism class.
// The class contains non-static methods, that have knowledge of the Automorphism instance on which they
// are invoked, and static methods, to which Automorphism instances may be passed as function arguments.
// Although we refrain from the final keyword, we treat all arrays as immutable. Thus, it is not a problem if
// different Automorphism instances refer to the same int[][] or if different int[][] refer to the same int[].
class Automorphism {

  // We encode the automorphism represented Automorphism instances by its values on the generators of F_4.
  // e.g. {{1,2,-1},{1,3},{1},{4}} encodes the automorphism "abA, ac, a, d, " of F_4 = <a, b, c, d> that maps
  // a to aba^{-1}, b to ac, c to a and d to d. The letters a, b, c, d correspond to the numbers 1, 2, 3, 4,
  // and negative numbers correspond to capital letters, i.e. inverse generators.
  int[][] values;

  // The elements of the diagonal subgroup (Z_2)^4 and the permutation subgroup S_4 of Aut(F_4), that
  // together generate the symmetry group G of the SDP problem, are simply represented by a String that
  // describes their action, e.g. "+-++" inverts the second generator of F_4, and "2134" permutes the
  // first two generators.
  // This is faster than treating them as Automorphism instances and less bloat than writing another class.
  // The order of the arrays must match the order in the accompanying MATLAB script.
  static String[] DIAGONALS = {"++++", "+++-", "++-+", "++--", "+-++", "+-+-", "+--+", "+---",
                               "-+++", "-++-", "-+-+", "-+--", "--++", "--+-", "---+", "----"};
  static String[] PERMUTATIONS = {"1234", "1243", "1324", "1342", "1423", "1432", "2134", "2143", "2314", "2341", "2413", "2431",
                                  "3124", "3142", "3214", "3241", "3412", "3421", "4123", "4132", "4213", "4231", "4312", "4321"};


  // Returns a new instance of the Automorphism class with the prescribed values field.
  public Automorphism(int[][] values) {
    this.values = values;
  }

  // Returns the identity Automorphism.
  static Automorphism identity() {
    int[][] values = new int[4][1];
    for (int i = 0; i < 4; i++) {
      values[i][0] = i+1; }
    return new Automorphism(values);
  }

  // Returns the Automorphism that sends the |i|-th generator of F_4 (its inverse if negative sign)
  // to the same element multiplied from right with the |j|-th generator (resp. its inverse). 1 <= i, j <= 4.
  static Automorphism elementaryAutomorphism(int i, int j) {
    int[][] values = new int[4][];
    for (int k = 0; k < 4; k++) {
      if (k == i-1) {
        values[k] = new int[] {k+1, j};
        continue; }
      if (k == -i-1) {
        values[k] = new int[] {-j, k+1};
        continue; }
      values[k] = new int[] {k+1};
    }
    return new Automorphism(values);
  }

  // Returns an array containing the elementary Automorphisms we use as generators.
  static Automorphism[] generators() {
    Automorphism[] result = new Automorphism[48];
    int pos = 0;
    for (int i = 1; i <= 4; i++) {
      for (int j = 1; j <= 4; j++) {
        if (i == j) continue;
        result[pos] = Automorphism.elementaryAutomorphism( i,  j); pos++;
        result[pos] = Automorphism.elementaryAutomorphism( i, -j); pos++;
        result[pos] = Automorphism.elementaryAutomorphism(-i,  j); pos++;
        result[pos] = Automorphism.elementaryAutomorphism(-i, -j); pos++;
      }
    }
    return result;
  }

  // Returns the product of the Automorphisms a and b
  static Automorphism product(Automorphism a, Automorphism b) {
    int[][] values = new int[4][];
    for (int i = 0; i < 4; i++) {
      values[i] = cook(a.values[i], b.values); }
    return new Automorphism(values);
  }

  // Helper method for product: Form a new word by concatenating the ingredient words
  // according to the recipe and taking care of letter cancellations
  private static int[] cook(int[] recipe, int[][] ingredients) {
    int[] temp = new int[100];   // we output first into this array of sufficiently large fixed size
    int pos = 0;                 // position in the output temp where we are currently working
    for (int i = 0; i < recipe.length; i++) {                    // for every step in the recipe
      int l = recipe[i];                                         // current letter, i.e., ingredient
      if (l > 0) {                                               // append l-th ingredient with cancellation
        for (int j = 0; j < ingredients[l-1].length; j++) {      // for each letter in ingredient
          if (pos > 0 && temp[pos-1]+ingredients[l-1][j] == 0) { // if it cancels preceding letter
            pos--; }                                             // delete that letter by moving back
          else {                                                 // if no cancellation, add letter and update position
            temp[pos] = ingredients[l-1][j];
            pos++; }
        }
      }
      else {                     // append inverse of l-th ingredient with cancellation in the same manner as above
        for (int j = ingredients[-l-1].length-1; j >= 0; j--) {
          if (pos > 0 && temp[pos-1]-ingredients[-l-1][j] == 0) {
            pos--; }
          else {
            temp[pos] = -ingredients[-l-1][j];
            pos++; }
        }
      }
    }
    int[] result = new int[pos]; // now copy the output into an array of correct length
    for (int i = 0; i < result.length; i++) {
      result[i] = temp[i]; }
    return result;
  }


  // Tests whether two instances of Automorphism encode the same element. Used by Hashtable.
  @Override
  public boolean equals(Object obj) {
    if (!(obj instanceof Automorphism)) return false;
    int[][] otherValues = ((Automorphism) obj).values;
    return Arrays.deepEquals(this.values, otherValues);
  }

  // Hashtable and HashSet use this method. Compatible with equals, i.e.,
  // Automorphism instances that encode the same group element give the same hash.
  @Override
  public int hashCode() {
    return Arrays.deepHashCode(this.values);
  }

  // Returns a (human readable) String representation of the Automorphism.
  // Different Automorphisms give different Strings, so the String can be used as a proxy.
  @Override
  public String toString() {
    String result = "";
    for (int i = 0; i < values.length; i++) {
      for (int j = 0; j < values[i].length; j++) {
        if (values[i][j] == 0) continue;
        if (values[i][j] < 0) {
          result += String.valueOf((char)(-values[i][j]-1 + 'A')); }
        if (values[i][j] > 0) {
          result += String.valueOf((char)(values[i][j]-1 + 'a')); }
      }
      result += ", ";
    }
    return result;
  }


  // The function isSmaller takes the values field of two Automorphisms and checks if the first Automorphism
  // comes before the second in our (arbitrary but fixed) total order: The generators and their inverses
  // are given the order of letterOrder; words in F_4 are ordered by length, then lexicographically;
  // Automorphisms are ordered by the values they assign to the generators.
  private static int[] letterOrder = {8, 6, 4, 2, 0, 1, 3, 5, 7}; // a < A < b < B < c < C < d < D
  static boolean isSmaller(Automorphism a, Automorphism b) {
    for (int i = 0; i < a.values.length; i++) {
      if (a.values[i].length < b.values[i].length) return false;
      if (a.values[i].length > b.values[i].length) return true;
      for (int j = 0; j < a.values[i].length; j++) {
        int comp = letterOrder[4+a.values[i][j]] - letterOrder[4+b.values[i][j]];
        if (comp < 0) return true;
        if (comp > 0) return false;
      }
    }
    return false;
  }


  // Returns the canonical representative for the point equivalence class of this automorphism.
  // The point equivalence class is the left coset with respect to the subgroup of Aut(F_4)
  // generated by the diagonal and the permutation subgroups. The canonical representative is
  // the minimal element in the equivalence class with respect to the order of isSmaller.
  // Instead of using isSmaller we calculate this element directly to increase speed:
  // Going through the the words values[i] word for word, letter by letter, we incrementally
  // build the representative as well as the symmetry T such that this*T = representative.
  Automorphism pointClassRepresentative() {
    int[][] representativeValues = new int[values.length][]; // the values field of the returned representative
    int[] translations = new int[4];                         // translation table representing the element T
    int[] takenLetters = new int[4];                         // generators (including inverses) already assigned in translations
    for (int i = 0; i < values.length; i++) {
      representativeValues[i] = new int[values[i].length];   // image of i-th generator, an element of F_4
      for (int j = 0; j < values[i].length; j++) {
        int letter = values[i][j];                           // the generator of F_4 occuring at position j in i-th image
        int sign = 1;                                        // -1 if values[i][j] is an inverse of a generator
        if (letter < 0) {
          letter = -letter;
          sign = -1; }
        if (translations[letter-1] == 0) {                   // if translations does not yet dictate a replacement
          for (int k = 1; k <= 4; k++) {                     // consider generators from best to worse
            if (takenLetters[k-1] == 0) {                    // take the first generator that is still available
              takenLetters[k-1] = 1;                         // mark it as unavailable for the future
              translations[letter-1] = k * sign;             // use it as replacement, with sign if that is better
              break;                                         // exit replacement-finding loop
            }
          }
        }
        representativeValues[i][j] = translations[letter-1]*sign; // write replacement into representative values
      }
    }   
    return new Automorphism(representativeValues);
  }


  // Returns the canonical representative for the distance equivalence class of the automorphism.
  // The distance equivalence class is the set of all automorphism obtained from a or its inverse
  // by left- and right multiplication with elements of the group generated by the diagonal and
  // the permutation subgroups. The canonical representative is the minimal element in this class
  // with respect to the order of isSmaller. To compute it we search for the the smallest representative
  // of the point equivalence classes into which the distance equivalence class decomposes.
  // This direct computation is much less memory intensive than a look-up in a hashtable.
  // The argument aInverse is allowed to differ from a^{-1} by left- or right multiplication with
  // an element in the symmetry group.
  static Automorphism distanceClassRepresentative(Automorphism a, Automorphism aInverse) {
    Automorphism currentBest = a;                    // the smallest representative found so far
    for (Automorphism e : a.diagonalActionOrbit()) { // diagonalOrbit computes slower, hence outer loop
      for (Automorphism f : e.permutationActionOrbit()) {
        Automorphism g = f.pointClassRepresentative();
        if (isSmaller(g, currentBest)) {
          currentBest = g; }
      }
    }
    for (Automorphism e : aInverse.diagonalActionOrbit()) {
      for (Automorphism f : e.permutationActionOrbit()) {
        Automorphism g = f.pointClassRepresentative();
        if (isSmaller(g, currentBest)) {
          currentBest = g; }
      }
    }
    return currentBest;
  }


  // Returns an array containing the automorphisms obtained from this one by left multiplication with
  // the elements of the diagonal subgroup, in order.
  Automorphism[] permutationActionOrbit() {
    Automorphism[] result = new Automorphism[PERMUTATIONS.length];
    for (int i = 0; i < PERMUTATIONS.length; i++) {
      int[][] permutedValues = new int[values.length][];
      for (int j = 0; j < values.length; j++) {
        permutedValues[j] = values[PERMUTATIONS[i].charAt(j)-'1']; }
      result[i] = new Automorphism(permutedValues);
    }
    return result;
  }

  // Returns an array containing the automorphisms obtained from this one by left multiplication with
  // the elements of the permutation subgroup, in order.
  // To save time we first compute the inverse of all generator images.
  Automorphism[] diagonalActionOrbit() {
    int[][] flippedValues = new int[values.length][]; // inverses of the generator images
    for (int i = 0; i < values.length; i++) {
      flippedValues[i] = new int[values[i].length];
      for (int j = 0; j < values[i].length; j++) {
        flippedValues[i][j] = -values[i][values[i].length-1-j];
      }
    }
    Automorphism[] result = new Automorphism[DIAGONALS.length];
    for (int i = 0; i < DIAGONALS.length; i++) {
      int[][] diagonalledValues = new int[values.length][];
      for (int j = 0; j < values.length; j++) {
        if (DIAGONALS[i].charAt(j) == '-') {
          diagonalledValues[j] = flippedValues[j]; }
        else {
          diagonalledValues[j] = values[j]; }
      }
      result[i] = new Automorphism(diagonalledValues);
    }
    return result;
  }

}
