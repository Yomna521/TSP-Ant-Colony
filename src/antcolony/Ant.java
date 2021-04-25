package antcolony;

import java.io.*;
import java.util.*;

public class Ant {

    private static int graphSize;
    private static double[][] costs = {    //Cities' distance settings
				{0, 13, 4, 6, 27, 1, 11, 3, 5, 6},
				{13, 0, 2, 2, 5, 6, 3, 4, 1, 1},
				{4, 2, 0, 4, 4, 2, 2, 9, 7, 9},
				{6, 2, 4, 0, 10, 11, 11, 14, 19, 17},
				{27, 5, 4, 10, 0, 22, 28, 18, 19, 33},
				{1, 6, 2, 11, 22, 0, 63, 15, 18, 10},
				{11, 3, 2, 11, 28, 63, 0, 17, 20, 20},
				{3, 4, 9, 14, 18, 15, 17, 0, 11, 10},
				{5, 1, 7, 19, 19, 18, 20, 11, 0, 10},
				{6, 1, 9, 17, 33, 10, 20, 10, 10, 0}
    };
    

    private static final Random rand = new Random(System.currentTimeMillis());

    private static double pheromoneLevels[][];
    

    /* Parameters */
    private static final int NO_OF_ANTS = 10;
    private static final int ITERATIONS = 150;
    private static final double ALPHA = 0.1;
    private static final double BETA = 1;
    private static final double RHO = 0.1;
    private static final double Q0 = 0.9;

    private static final ArrayList<ArrayList<Integer>> ants = new ArrayList<>(NO_OF_ANTS);
    private static ArrayList<Integer> globalBestAnt = new ArrayList<>();
    private static double globalBestTourLength;
    //private static double iterationBestTourLength;

    public static void main(String[] args) throws FileNotFoundException {
        costs = DataParser.parse("C:\\Users\\Yumna\\Downloads\\TSP-Search-Algorithm-master\\TSP-Search-Algorithm-master\\DATA\\Cincinnati.tsp");
        graphSize = costs.length;
        System.out.println(Arrays.deepToString(costs));
        long startTime = System.currentTimeMillis();
        ArrayList<Integer> bestTour = aco(costs);
        long endTime   = System.currentTimeMillis();
        long totalTime = endTime - startTime;
        System.out.println("------------------------ Best Tour ----------------------");
        System.out.println(bestTour);
        System.out.println("Cost: " + computeTourLength(globalBestAnt));
        System.out.println("Total time to execute: "+totalTime+" ms");
    }

    /**
    Ant Colony System for TSP. Implements TSP tours as a list of cities. Each tour is
    logically closed at the point of computing its cost.
     **/
    @SuppressWarnings("empty-statement")
    private static ArrayList<Integer> aco(double graph[][]) {

        /* Initialization */

        double initialPheromoneLevel = 1 / (graphSize * nearestNeighbor(graph));
        pheromoneLevels = new double[graph.length][graph.length];

        
        for (int i = 0; i < NO_OF_ANTS; i++) {
            
            ants.add(new ArrayList<Integer>());
        }


        //System.out.println(Arrays.deepToString(pheromoneLevels));
        /* Search - Using ACS */
        ArrayList<Integer> iterationBestAnt;
        // Pheromone level
        for (int i = 0; i < graphSize; i++) {  // using graphSize for quality assurance :|
            for (int j = 0; j < graphSize; j++) {
                if (i == j) {
                    pheromoneLevels[i][j] = Double.NaN;
                } else {
                    pheromoneLevels[i][j] = initialPheromoneLevel;
                }
            }
        }
        for (int iter = 0; iter < ITERATIONS; iter++) {

            // Initialize Ants
            ArrayList<Integer> startPoints = new ArrayList<>();
            for (int i = 0; i < NO_OF_ANTS; i++) {
                // Each ant starts at a random node, and at most one node on each node
                int startPoint;
                while (startPoints.contains(startPoint = rand.nextInt(graphSize))) ;
                startPoints.add(startPoint);
                ants.get(i).add(startPoint);
            }
            for (int j = 1; j < graphSize; j++) {
                int r = j - 1;
                for (int k = 0; k < NO_OF_ANTS; k++) {

                    /* Exploitation */
                    double maxChoice = Double.MIN_VALUE;
                    int maxIndex = 0;
                    for (int u = 0; u < graphSize; u++) {
                        // Implicit TabuList - avoid adding nodes which are already part of the tour
                        if (ants.get(k).contains(u)) continue;

                        double desirable = pheromoneLevels[ants.get(k).get(r)][u]
                                * Math.pow(1 / graph[ants.get(k).get(r)][u], BETA);
                        if (desirable > maxChoice) {
                            maxChoice = desirable;
                            maxIndex = u;
                        }
                    }
                    double q = rand.nextDouble();
                    if (q < Q0) {
                        if(ants.get(k).contains(maxIndex)){
                            throw new ArrayIndexOutOfBoundsException(ants.get(k).toString() + ": " + maxIndex);
                        }
                        ants.get(k).add(maxIndex);

                    } else {
                        /* Biased exploration */
                        double totalBias = 0;
                        for (int u = 0; u < graphSize; u++) {
                            // Implicit TabuList - avoid adding nodes which are already part of the tour
                            if (ants.get(k).contains(u)) continue;
                            totalBias += pheromoneLevels[ants.get(k).get(r)][u]
                                    * Math.pow(1 / graph[ants.get(k).get(r)][u], BETA);
                        }

                        maxChoice = Double.MIN_VALUE;
                        maxIndex = -1;
                        for (int s = 0; s < graphSize; ++s) {
                            // Implicit TabuList - avoid adding nodes which are already part of the tour
                            if (ants.get(k).contains(s)){
                                continue;
                            }

                            double numerator = pheromoneLevels[ants.get(k).get(r)][s]
                                    * Math.pow(1 / graph[ants.get(k).get(r)][s], BETA);

                            double probabilityK = numerator / totalBias;
                            //System.out.println("r, s: " + r + " " + s + "probability: " + probabilityK);
                            if (probabilityK > maxChoice) {
                                maxChoice = probabilityK;
                                maxIndex = s;
                            }
                        }
                        if(ants.get(k).contains(maxIndex)){
                            throw new ArrayIndexOutOfBoundsException("Duplicate!\n" + ants.get(k).toString() + ": " + maxIndex);
                        }
                        ants.get(k).add(maxIndex);

                        // Local update
                        pheromoneLevels[ants.get(k).get(j-1)][maxIndex] = (1-RHO) * pheromoneLevels[ants.get(k).get(j-1)][maxIndex]
                                + RHO * initialPheromoneLevel;
                    }
                } // Ants
            }
            iterationBestAnt = getIterationBestAnt();
            double iterationBestAntCost = computeTourLength(iterationBestAnt);

            if (globalBestAnt.isEmpty()){
                globalBestAnt = new ArrayList<>(getIterationBestAnt());
                globalBestTourLength =  computeTourLength(globalBestAnt);
                System.out.println(globalBestAnt + " " + globalBestTourLength);
                globalUpdate(globalBestAnt);
            } else if(iterationBestAntCost < globalBestTourLength){
                globalBestAnt = new ArrayList<>(iterationBestAnt);
                globalBestTourLength =  computeTourLength(globalBestAnt);
                System.out.println(globalBestAnt + " " + globalBestTourLength);
                globalUpdate(globalBestAnt);
            }

            globalBestAnt = twoOptLocalSearch(globalBestAnt);

            // Clear ants information
            for(ArrayList<Integer> ant: ants){
                ant.clear();
            }
        } // Iterations

        return globalBestAnt;
    }

    public static ArrayList<Integer> twoOptMove(ArrayList<Integer> tour, int j, int k){
        ArrayList<Integer> newTour = new ArrayList<>(tour.subList(0, j));
        ArrayList<Integer> reverseTour = new ArrayList<>(tour.subList(j, k));
        Collections.reverse(reverseTour);
        newTour.addAll(reverseTour);
        ArrayList<Integer> tourEnd = new ArrayList<>(tour.subList(k, tour.size()));
        newTour.addAll(tourEnd);

        return newTour;
    }


    public static ArrayList<Integer> twoOptRandom(ArrayList<Integer> tour){
        int nCities = tour.size();
        int j = rand.nextInt(nCities);
        int k = rand.nextInt(nCities);
        while (j == k) k = rand.nextInt(nCities);

        // Ensure that j is less than k. That is, if j > k, swap j and k
        int tempCut = j > k ? j : k;
        if (tempCut == j) {
            j = k;
            k = tempCut;
        }
        return twoOptMove(tour, j, k);
    }

    public static ArrayList<Integer> twoOptLocalSearch(ArrayList<Integer> tour){
        for (int i = 0; i < tour.size(); i++) {
            for (int j = 0; j < tour.size(); j++) {
                if(i < j){
                    ArrayList<Integer> newTour = twoOptMove(tour, i, j);
                    if (computeTourLength(newTour) < computeTourLength(tour)){
                        return newTour;
                    }
                }
            }
        }
        return tour;
    }

    /**
     * Only occurs for the best tour
     * @param tour the global best tour corresponding to the best ant
     *
     */
    public static void globalUpdate(ArrayList<Integer> tour){
        double tourLength = computeTourLength(tour);
        for (int i = 0; i < tour.size()-1; i++) {
            pheromoneLevels[tour.get(i)][tour.get(i+1)] = (1-ALPHA) * pheromoneLevels[tour.get(i)][tour.get(i+1)]
                    + ALPHA * (1/tourLength);
        }

        // Close tour
        pheromoneLevels[tour.get(tour.size()-1)][tour.get(0)] = (1-ALPHA) * pheromoneLevels[tour.get(tour.size()-1)][tour.get(0)]
                + ALPHA * (1/tourLength);
    }



    private static ArrayList<Integer> getIterationBestAnt(){
        double minTourLength = Double.MAX_VALUE;
        int minIndex = 0;
        int i = 0;
        for(ArrayList<Integer> ant: ants){
            double tourLength = computeTourLength(ant);
            if (tourLength < minTourLength){
                minIndex = i;
            }
            ++i;
        }
        return ants.get(minIndex);
    }


    public static double nearestNeighbor(double graph[][]) {
        Random ran = new Random(System.currentTimeMillis());
        int start = ran.nextInt(graph.length);
        ArrayList<Integer> tour = new ArrayList<>();
        tour.add(start);
        double optimalTourLength = 0;
        int i = 0;
        int minNeighborIndex = 0;
        while (tour.size() < graph.length) {
            double minNeighborLength = Double.MAX_VALUE;
            for (int j = 0; j < graph.length; j++) {
                if (graph[i][j] < minNeighborLength && !tour.contains(j)) {
                    minNeighborLength = graph[i][j];
                    minNeighborIndex = j;
                }
            }
            tour.add(minNeighborIndex);
            optimalTourLength += minNeighborLength;
            i = minNeighborIndex;
        }

        optimalTourLength += graph[minNeighborIndex][0]; // Close the tour
        //System.out.println(tour);
        return optimalTourLength;
    }

    public static double computeTourLength(ArrayList<Integer> tour){
        double tourLength = 0;
        for (int i = 0; i < tour.size()-1; i++) {
            tourLength += costs[tour.get(i)][tour.get(i+1)];
        }
        tourLength += costs[tour.get(tour.size()-1)][tour.get(0)]; // Close tour;
        return tourLength;
    }
    
}