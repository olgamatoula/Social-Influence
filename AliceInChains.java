package gr.james.socialinfluence.tournament.players;

import gr.james.socialinfluence.game.Move;
import gr.james.socialinfluence.game.players.Player;
import gr.james.socialinfluence.graph.Edge;
import gr.james.socialinfluence.graph.Graph;
import gr.james.socialinfluence.graph.Vertex;
import gr.james.socialinfluence.graph.algorithms.Degree;
import gr.james.socialinfluence.graph.algorithms.Dijkstra;
import gr.james.socialinfluence.graph.algorithms.PageRank;
import gr.james.socialinfluence.graph.collections.GraphState;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Set;

public class AliceInChains extends Player {

	HashMap<Vertex, Double> sums;
	String graphName = "";
	int initialClique = 0;
	int stepEdges = 0;
	int numOfVertices = 0;
	int numOfClusters = 0;

	@Override
	public void getMove() {
		
		long startTime = System.currentTimeMillis();
		String graphType = this.g.getMeta();
		boolean verbose = !this.d.getTournament();
		//num of moves o arithmos twn epitreptwn kinhsewn
		int numOfMoves = this.d.getNumOfMoves();
		// numOfVertices panta monos arithmos
		int numOfNodes = this.g.getVerticesCount();
		Move m = new Move();
		parseGraphType(graphType);
		
		if (numOfMoves < numOfNodes) {
			if (graphName.equals("TwoWheels")) {
				// to seed sto generate()
				int n = (numOfNodes + 1) / 2;

				if (verbose) {
					System.out.println("GraphType: " + graphType
							+ " with initial seed: " + n + "\n"
							+ "Number of moves: " + numOfMoves
							+ " Number of vertices: " + numOfVertices);
				}

				GraphState degrees = Degree.execute(g, true);

				int[] largestDegrees = new int[3];
				Vertex[] largestDegreeVertices = new Vertex[3];
				for (int i = 0; i < 3; i++) {
					largestDegrees[i] = 0;
					largestDegreeVertices[i] = this.g.getVertexFromId(2);
				}

				//find the 3 largest values and their vertex position
				for (Vertex v : degrees.keySet()) {
					Double val = new Double(degrees.get(v).doubleValue());
					if (val.intValue() > largestDegrees[0]) {
						largestDegrees[2] = largestDegrees[1];
						largestDegrees[1] = largestDegrees[0];
						largestDegrees[0] = val.intValue();
						largestDegreeVertices[2] = largestDegreeVertices[1];
						largestDegreeVertices[1] = largestDegreeVertices[0];
						largestDegreeVertices[0] = v;
					} else if (val.intValue() > largestDegrees[1]) {
						largestDegrees[2] = largestDegrees[1];
						largestDegrees[1] = val.intValue();
						largestDegreeVertices[2] = largestDegreeVertices[1];
						largestDegreeVertices[1] = v;
					} else if (val.intValue() > largestDegrees[2]) {
						largestDegrees[2] = val.intValue();
						largestDegreeVertices[2] = v;
					}
				}

				if (verbose) {
					System.out.println();
					System.out
							.println("The 3 largest nodes according to degree centrality. ");
					for (int i = 0; i < 3; i++) {
						System.out.println("Node " + largestDegreeVertices[i]
								+ " with Degree " + largestDegrees[i]);
					}
				}

				Vertex cg = null, cw1 = null, cw2 = null;
				if (n <= 4) {
					cg = largestDegreeVertices[0];
					cw1 = largestDegreeVertices[1];
					for (Vertex v : degrees.keySet()) {
						if (!(v == cg || v == cw1)) {
							if (!hasEdge(v, cw1)) {
								cw2 = v;
								break;
							}
						}
					}
				} else if (n == 5 && ((numOfMoves == 3) || (numOfMoves == 4))) {
					cg = largestDegreeVertices[0];
					List<Vertex> centers = new ArrayList<Vertex>();
					centers.add(cw1);
					centers.add(cw2);
					for (Vertex v : this.g.getVertices()) {
						if ((!hasEdge(v, cg)) && (v != cg)) {
							cw1 = v;
							break;
						}
					}
					for (Vertex v : this.g.getVertices()) {
						if ((!hasEdge(v, cg)) && (v != cg) && (v != cw1)) {
							cw2 = v;
							break;
						}
					}
				} else if (n < 7) {
					cg = largestDegreeVertices[0];
					cw1 = largestDegreeVertices[1];
					cw2 = largestDegreeVertices[2];

				} else if (n > 7) {
					cg = largestDegreeVertices[2];
					cw1 = largestDegreeVertices[0];
					cw2 = largestDegreeVertices[1];
				} else {
					Set<Edge> edges = largestDegreeVertices[0].getOutEdges();
					boolean found1 = false, found2 = false;
					for (Edge e : edges) {
						if (e.getTarget() == largestDegreeVertices[1])
							found1 = true;
						else if (e.getTarget() == largestDegreeVertices[2])
							found2 = true;
						if (found1 && found2)
							break;
					}
					if (found1 && found2) {
						cg = largestDegreeVertices[0];
						cw1 = largestDegreeVertices[1];
						cw2 = largestDegreeVertices[2];
					} else if (found1) {
						cg = largestDegreeVertices[1];
						cw1 = largestDegreeVertices[0];
						cw2 = largestDegreeVertices[2];
					} else if (found2) {
						cg = largestDegreeVertices[2];
						cw1 = largestDegreeVertices[0];
						cw2 = largestDegreeVertices[1];
					}
				}

				double cgBudget = (double) 1;
				double cwBudget = ((double) (n - 1) / 3.883495146) * cgBudget;

				switch (numOfMoves) {
				case 1:
					m.putVertex(cg, 1);
					break;
				case 2:
					if (n <= 4) {
						m.putVertex(cg, 0.75);
						m.putVertex(cw1, 0.25);
					} else {
						m.putVertex(cw1, 0.5);
						m.putVertex(cw2, 0.5);
					}
					break;
				case 3:
					if (n == 6) {
						m.putVertex(cg, 1);
						m.putVertex(cw1, 1.2);
						m.putVertex(cw2, 1.2);
					} else if (n == 5) {
						m.putVertex(cg, 1.7);
						m.putVertex(cw1, 1);
						m.putVertex(cw2, 1);
					} else {
						m.putVertex(cg, cgBudget);
						m.putVertex(cw1, cwBudget);
						m.putVertex(cw2, cwBudget);
					}
					break;
				default:
					/*
					 * movesLeft is the number of moves that need to be defined for each half circle
					 * after the first 3 moves are defined
					 */
					int nodesOnEachSide = n - 2;
					int movesLeft = numOfMoves - 3;
					int movesLeft1 = 0,
					movesLeft2 = 0;
					int stepsTNP1 = 0,
					stepsTNP2 = 0;
					if (movesLeft % 2 == 0) {
						movesLeft1 = movesLeft / 2;
						movesLeft2 = movesLeft1;
						stepsTNP1 = (int) Math.round((double) nodesOnEachSide
								/ (movesLeft1 + 1));
						if (stepsTNP1 == 0)
							stepsTNP1 = 1;
						stepsTNP2 = stepsTNP1;
					} else {
						movesLeft1 = (movesLeft + 1) / 2;
						movesLeft2 = movesLeft1 - 1;
						stepsTNP1 = (int) Math.round((double) nodesOnEachSide
								/ (movesLeft1 + 1));
						if (stepsTNP1 == 0)
							stepsTNP1 = 1;
						stepsTNP2 = (int) Math.round((double) nodesOnEachSide
								/ (movesLeft2 + 1));
					}
					if (movesLeft1 > ((double) nodesOnEachSide / 2)) {
						stepsTNP1 = 2;
						if (movesLeft2 > ((double) nodesOnEachSide / 2)) {
							stepsTNP2 = 2;
						}
					}

					double nBudget = (double) 1;
					double ratio = ((double) numOfMoves / numOfVertices);
					//	            double ratio = 1;
					cgBudget = 1.27 * nBudget;
					//	            cwBudget = ((double) (n - 1) / 3.883495146 ) * cgBudget;
					//	            if((ratio > (double) 0.75) && n != 4){
					cwBudget = ((double) ((n - 1) / (ratio * 10))) * cgBudget;
					//	            }

					m.putVertex(cg, cgBudget);
					m.putVertex(cw1, cwBudget);
					m.putVertex(cw2, cwBudget);

					Vertex source = cg;
					Vertex previous = cw2;
					Vertex ppt = null;
					for (int i = 0; i < movesLeft1; i++) {
						for (int j = 0; j < stepsTNP1; j++) {
							ppt = chooseVertex(source, previous, cw1);
							previous = source;
							source = ppt;
						}
						m.putVertex(ppt, nBudget);
					}
					source = cg;
					previous = cw1;
					ppt = null;
					for (int i = 0; i < movesLeft2; i++) {
						for (int j = 0; j < stepsTNP2; j++) {
							ppt = chooseVertex(source, previous, cw2);
							previous = source;
							source = ppt;
						}
						m.putVertex(ppt, nBudget);
					}
					int sizeOfMove = m.getVerticesCount();
					//an perissepsan kinhseis, pare prwta ayta poy einai pio konta sto kentro
					//kai meta dialekse random nodes
					//		        if(sizeOfMove<numOfMoves){
					//		        	Set<Edge> cgEdges = cg.getOutEdges();
					//		        	for (Edge e:cgEdges){
					//		        		if(!m.containsVertex(e.getTarget())&&(sizeOfMove<numOfMoves)){
					//		        			m.putVertex(e.getTarget(), nBudget);
					//		        			sizeOfMove++;
					//		        		}
					//		        	}
					//		        }
					boolean side = true;
					while (sizeOfMove < numOfMoves) {
						Vertex v = this.g.getRandomVertex();
						if (!m.containsVertex(v)) {
							if (side) {
								if (hasEdge(v, cw1)) {
									m.putVertex(v, nBudget);
									sizeOfMove++;
									side = !side;
								}
							} else {
								if (hasEdge(v, cw2)) {
									m.putVertex(v, nBudget);
									sizeOfMove++;
									side = !side;
								}
							}
						}
					}
					break;
				}
			} // TwoWheels ends here
			else if (graphName.equals("BarabasiAlbert")) {			
				GraphState pagerank = PageRank.execute(g, 0.0);
				HashMap<Vertex, Double> eccentricities = calculateEccentricity(this.g);
				Double minPathSums = findMin(sums,this.g);
				
				Double[][] sortedPagerank =  sortToNumOfMoves(pagerank, numOfNodes, false); // descending
				List<Double> PagerankID = new ArrayList<Double>(Arrays.asList(sortedPagerank[0]));
				Double[][] sortedEc = sortToNumOfMoves(eccentricities, numOfNodes, true); //ascending
				Double[][] sortedSums = sortToNumOfMoves(sums, numOfMoves, true); //ascending
				
				
				double weight = 1.0;
				double ratio = (double) numOfNodes/numOfMoves;
				Vertex v0,v1;
				
				int i = 0;
				int sizeOfMove = 0;
				boolean EdgeinMove = false;
				boolean EdgeinMove2 = false;
				
				if (numOfMoves == 1) {
					v0 = this.g.getVertexFromId((int) ((double) PagerankID.get(0)));
					v1 = this.g.getVertexFromId((int) ((double) sortedEc[0][0]));
					if ((Double.compare(eccentricities.get(v0), sortedEc[1][0]) == 0) || (0.8 * pagerank.get(v0) > pagerank.get(v1))){
						m.putVertex(v0, 1);
					} else {
						m.putVertex(v1, 1);
					}
				} else { 
				
					while (sizeOfMove < numOfMoves) {
		
						i = numOfMoves - sizeOfMove;
						EdgeinMove = false;
						EdgeinMove2 = false;
						v0 = this.g.getVertexFromId((int) ((double) PagerankID.get(0)));		
						v1 = this.g.getVertexFromId((int) ((double) PagerankID.get(1)));
						weight = (pagerank.get(v0));
						
						if ((pagerank.get(v0) - pagerank.get(v1)) < (0.1 * (pagerank.get(v0)))) {
							
							Set<Edge> edges = v0.getOutEdges();
							for (Edge e: edges) {
								if (m.containsVertex(e.getTarget())) {
									EdgeinMove = true;
									break;
								} 						
							}
							if (EdgeinMove) {
								Set<Edge> edges2 = v1.getOutEdges();
								for (Edge e : edges2) {
									if (m.containsVertex(e.getTarget())) {
										EdgeinMove2 = true;
										break;
									}
								}
							}
							if ( eccentricities.get(v0) > eccentricities.get(v1) ) {
								m.putVertex(v1, weight);
								PagerankID.remove(1);
							}else if ((!EdgeinMove) || (EdgeinMove2)) {
								weight = weight * ((i / 100) + 1);
								m.putVertex(v0, weight);
								PagerankID.remove(0);
							} else if (!EdgeinMove2) {
								weight = weight * ((i / 100) + 1);
								m.putVertex(v1, weight);
								PagerankID.remove(1);
							}				
							
						} else if ((pagerank.get(v0) - pagerank.get(v1)) < (0.2 * (pagerank.get(v0)))) {
							weight = weight * ((i / 100) + 1);
							m.putVertex(v0, weight);
							PagerankID.remove(0);
						
						} else if ((pagerank.get(v0) - pagerank.get(v1)) > ( 0.3 * (pagerank.get(v0)))) {
							weight = weight / ((i / 100) + 1);
							m.putVertex(v0, weight);
							PagerankID.remove(0);
						} else {
							m.putVertex(v0, weight);
							PagerankID.remove(0);
						}
						
						sizeOfMove = m.getVerticesCount();
					}	
				}
			} // BarabasiAlbert ends here
			else if (graphName.equals("BarabasiAlbertCluster")) {
				GraphState pagerank = PageRank.execute(g, 0.0);

				HashMap<Vertex, Double> eccentricities = calculateEccentricity(this.g);
				//after calculateEccentricity, sums have been calculated as well
				Double minPathSums = findMin(sums, this.g);

				Double[][] sortedPagerank = sortToNumOfMoves(pagerank,
						numOfNodes, false); // descending
				Double[][] sortedSums = sortToNumOfMoves(sums, numOfNodes, true); //ascending
				double ratio = (double) numOfMoves/numOfNodes;
				int p = 0;
				double weight = 1.0;
				//just put the top pagerank values and vertices in case we run out of time;
				while (m.getVerticesCount() < numOfMoves) {
					if(ratio <=  0.5)
						weight = sortedPagerank[1][p];
					else
						weight = Math.log10(10 * sortedPagerank[1][p]); 
					m.putVertex(
							this.g.getVertexFromId((int) (double) sortedPagerank[0][p]),
							weight);
					p++;
				}
				this.movePtr.set(m);

				Vertex[][] clusterIds = new Vertex[numOfClusters][numOfVertices];
				for (int i = 0; i < numOfClusters; i++) {
					clusterIds[i][0] = this.g
							.getVertexFromId((int) (double) sortedSums[0][i]);
				}

				for (int i = 0; i < numOfClusters; i++) {
					for (int j = 0; j < numOfVertices; j++) {
						Vertex next = clusterIds[i][j];
						Set<Edge> edges = next.getOutEdges();
						for (Edge e : edges) {
							Vertex target = e.getTarget();
							if (Double.compare(sums.get(target), minPathSums) == 0)
								continue;
							boolean exists = false;
							int k = 0;
							while (k < numOfVertices) {
								k++;
								if (clusterIds[i][k] == null) {
									break;
								} else if (target == clusterIds[i][k]) {
									exists = true;
									break;
								}
							}
							if (!exists) {
								clusterIds[i][k] = target;
							}
						}
					}
				}

				Double[][][] sortedClusters = new Double[numOfClusters][2][numOfVertices];
				int movesPerCluster = (int) numOfMoves / numOfClusters;
				int remainingMoves = numOfMoves - movesPerCluster
						* numOfClusters;
				for (int i = 0; i < numOfClusters; i++) {
					HashMap<Vertex, Double> tmps = new HashMap<Vertex, Double>();
					for (int j = 0; j < numOfVertices; j++) {
						Vertex temp = clusterIds[i][j];
						tmps.put(temp, pagerank.get(temp));
					}
					sortedClusters[i] = sortToNumOfMoves(tmps, numOfVertices,
							false);
				}

				m.clear();

				for (int i = 0; i < numOfClusters; i++) {
					for (int j = 0; j < movesPerCluster; j++) {
						if(ratio <=  0.5)
							weight = sortedClusters[i][1][j];
						else
							weight = Math.log10(10 * sortedClusters[i][1][j]); 
						m.putVertex(
								this.g.getVertexFromId((int) (double) sortedClusters[i][0][j]),
								weight);
					}
				}
				while (remainingMoves > 0) {
					Vertex v = this.g.getVertexFromId((int) (double) sortedClusters[(numOfClusters - remainingMoves)][0][movesPerCluster]);
					if(ratio <=  0.5)
						weight = sortedClusters[(numOfClusters - remainingMoves)][1][movesPerCluster];
					else
						weight = Math.log10(10 * sortedClusters[(numOfClusters - remainingMoves)][1][movesPerCluster]);
					m.putVertex(v,weight);
					remainingMoves--;
				}
			} // BarabasiAlbertCluster ends here
			else {
				GraphState pagerank = PageRank.execute(g, 0.85);

				Double[][] sortedPagerank = sortToNumOfMoves(pagerank,
						numOfNodes, false); // descending

				double weight = 1.0;

				Vertex v0;
				for (int i = 0; i < numOfMoves; i++) {
					v0 = this.g
							.getVertexFromId((int) ((double) sortedPagerank[0][i]));
					//				weight = pagerank.get(v0);
					weight = Math.log10(10 * pagerank.get(v0));
					m.putVertex(v0, weight);
				}
			} //random graph stuff end here
		} else if (numOfMoves == numOfNodes){
			for (Vertex v : this.g.getVertices()) {
				m.putVertex(v, 1);
			}
		}
		if (verbose){
			System.out.println("EXECUTION TIME(ms): " + ((long) System.currentTimeMillis()-startTime));
		}
		this.movePtr.set(m);
	} //getMove ends here
	private Vertex chooseVertex(Vertex source, Vertex previous, Vertex circleCenter) {
		Set<Edge> edges = source.getOutEdges();
		for (Edge e: edges) {
			if ((e.getTarget() != previous) && (e.getTarget() != circleCenter) && (hasEdge(e.getTarget(),circleCenter)))
				return e.getTarget();
		}
		return null;
	}
	
	private boolean hasEdge(Vertex source, Vertex target) {
		for (Edge e: source.getOutEdges()) {
			if (e.getTarget() == target)
				return true;
		}
		return false;
	}
	
	public void parseGraphType(String graphType){
		String[] tokens = graphType.split(",");
		int size = tokens.length;
//		for (int i = 0; i < size; i++){
//			System.out.println(tokens[i]);
//		}
		if(size > 0){
			graphName = tokens[0];			
			if (graphName.equals("BarabasiAlbert")||graphName.equals("BarabasiAlbertCluster")
					||graphName.equals("Path")||graphName.equals("TwoWheels")) {
				if (size == 2) { //TwoWheels
					String[] twTokens = tokens[1].split("=");
					numOfVertices = Integer.parseInt(twTokens[1]);
				} else if (size == 3) { // Path
					String[] pathTokens = tokens[1].split("=");
					numOfVertices = Integer.parseInt(pathTokens[1]);
				} else if (size == 6) { //BarabasiAlbert
					String[] baTokens = tokens[1].split("=");
					numOfVertices = Integer.parseInt(baTokens[1]);
					baTokens = tokens[2].split("=");
					initialClique = Integer.parseInt(baTokens[1]);
					baTokens = tokens[3].split("=");
					stepEdges = Integer.parseInt(baTokens[1]);
					//bugfix for cluster
					if (tokens[4].length()>3){ //e.g. a=1.000000
						baTokens = tokens[5].split("=");
						numOfClusters = Integer.parseInt(baTokens[1]);
					}
					
				} else if (size == 7) {
					String[] bacTokens = tokens[1].split("=");
					numOfVertices = Integer.parseInt(bacTokens[1]);
					bacTokens = tokens[2].split("=");
					initialClique = Integer.parseInt(bacTokens[1]);
					bacTokens = tokens[3].split("=");
					stepEdges = Integer.parseInt(bacTokens[1]);
					bacTokens = tokens[6].split("=");
					numOfClusters = Integer.parseInt(bacTokens[1]);
				}
			}
		}		
	}

	public HashMap<Vertex,Double> calculateEccentricity(Graph g){
		sums = new HashMap<Vertex,Double>();
		HashMap<Vertex, Double> eccentricities = new HashMap<Vertex, Double>();
		for (Vertex v : g.getVertices()){
			HashMap<Vertex,Double> hm = Dijkstra.execute(g, v);		
			for (Vertex u : g.getVertices() ){
				double val = 0;
				if (sums.get(v) != null) {
					val = sums.get(v) + hm.get(u);
				} else 
					val = hm.get(u);
				sums.put(v, val);
			}
			// ypologizei to eccentricity gia kathe komvo
			eccentricities.put(v, findMax(hm, g));
		}
		return eccentricities;
	}
	/*
	 * returns the maximum distance from a node to any other node
	 */
	public static Double findMax(HashMap<Vertex,Double> vectorDistances, Graph g){
		double max = -1;
		for (Vertex v : g.getVertices()){
			if (vectorDistances.get(v) > max) {
				max = vectorDistances.get(v);
			}
		}
		return max;
	}
	/*
	 * returns the minimum sum of path steps from any node to all the others 
	 */
	public static Double findMin(HashMap<Vertex, Double> sums, Graph g){
		double min = Double.POSITIVE_INFINITY;
		
		for (Vertex v : g.getVertices()){
			if (sums.get(v) < min) {
				min = sums.get(v);
			}
		}
		return min;
	}
	
	/*
	 * sort by type
	 * if type is true ascending, else descending
	 */
	public static Double[][] sortToNumOfMoves(HashMap<Vertex, Double> ec,int size, boolean type){
//		int size = ec.size();
		Double[] sortedVertices = new Double[size];
		Double[] sortedDoubles = new Double[size];
		
     	if (type){
     		for (int i = 0; i < size; i++)
     			sortedDoubles[i] = Double.POSITIVE_INFINITY;
     	} else {
     		for (int i = 0; i < size; i++)
     			sortedDoubles[i] = -1D;
     	}
     	
		for (Vertex v : ec.keySet()) {
	       	Double val = new Double(ec.get(v).doubleValue());
	        
       	for (int i = 0; i < size; i++){
       		if (type) {
				if (val < sortedDoubles[i]) {
					for (int j = size - 1; j > i; j--) {
						sortedDoubles[j] = sortedDoubles[j - 1];
						sortedVertices[j] = sortedVertices[j - 1];
					}
					sortedDoubles[i] = val;
					sortedVertices[i] = (double) v.getId();
					break;
				}
			} else {
				if (val > sortedDoubles[i]) {
					for (int j = size - 1; j > i; j--) {
						sortedDoubles[j] = sortedDoubles[j - 1];
						sortedVertices[j] = sortedVertices[j - 1];
					}
					sortedDoubles[i] = val;
					sortedVertices[i] = (double) v.getId();
					break;
				}
			}
       	}
		}
		Double[][] multi = new Double[2][];
		
		multi[0] = sortedVertices;
		multi[1] = sortedDoubles;
		return multi;
	}


}
