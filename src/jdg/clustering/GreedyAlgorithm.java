package jdg.clustering;
import jdg.graph.AdjacencyListGraph;
import jdg.graph.Node;

/**
 * This class provides an implementation of the Greedy algorithm for Community detection.
 * 
 * @author Luca Castelli Aleardi (INF421, 2017)
 */
public class GreedyAlgorithm extends CommunityDetection {
	
	public AdjacencyListGraph graph;
	public int nEdges;
	public int[][] matrixB;
	public int[] arrayK;
	public int[] communities;

	/**
	 * Initialize the parameters of the Louvain's algorithm
	 */
	public GreedyAlgorithm() {}
	
	
	/**
	 * This class lets us save during the loop the parameters which increase the modularity.
	 * At the end of the loop we know that we have to merge the two communities 'firstCommunityMaxModularityIncrease'
	 * and 'secondCommunityMaxModularityIncrease' to improve modularity.
	 * @author ziyed
	 *
	 */
	public class MaxModularityIncreaseParams{
		public int firstCommunityMaxModularityIncrease;
		public int secondCommunityMaxModularityIncrease;
		public double maxModularityIncrease;
		
		public MaxModularityIncreaseParams() {
			this.firstCommunityMaxModularityIncrease = 0;
			this.secondCommunityMaxModularityIncrease = 0;
			this.maxModularityIncrease = -2;
			
		}
	}
	
	
	/**
	 * This method finds the 'MaxModularityIncreaseParam' which maximizes the increase of modularity by merging two communities.
	 * @return maxModularityIncreaseParams
	 */
	public MaxModularityIncreaseParams findMaxModularityIncreaseParams() {
		int nodeCommunity;
		int neighbourCommunity;
		
		double modularityIncrease = 0;
		MaxModularityIncreaseParams maxModularityIncreaseParams  = new MaxModularityIncreaseParams();
		
		for(Node node: this.graph.vertices) {
			nodeCommunity = this.communities[node.index];
			for(Node neighbour: node.neighbors) {
				neighbourCommunity = this.communities[neighbour.index];
				if (nodeCommunity != neighbourCommunity) {
					modularityIncrease = (this.matrixB[nodeCommunity][neighbourCommunity] - this.arrayK[nodeCommunity] * this.arrayK[neighbourCommunity]/(2. * this.nEdges))/(this.nEdges);
					if (modularityIncrease > maxModularityIncreaseParams.maxModularityIncrease) {
						maxModularityIncreaseParams.maxModularityIncrease = modularityIncrease;
						maxModularityIncreaseParams.firstCommunityMaxModularityIncrease = nodeCommunity;
						maxModularityIncreaseParams.secondCommunityMaxModularityIncrease = neighbourCommunity;
					}
				}
			}
		}
		
		return maxModularityIncreaseParams;
	}
	
	
	
	/**
	 * This method returns a partition of a network of size 'n' into communities,
	 * computed by the Louvain's algorithm <p>
	 * <p>
	 * Remarks:<p>
	 * -) the nodes of the networks are numbered 0..n-1
	 * -) the graph is partitioned into 'k' communities, which have numbers 0..k-1
	 * 
	 * @param graph  the input network (adjacency list representation)
	 * @return an array of size 'n' storing, for each vertex, the index of its community (a value between 0..,k-1)
	 */
	public int[] computeClusters(AdjacencyListGraph graph) {
		long startTime = System.nanoTime();
		long endTime;
		this.graph = graph;
		int nVertices = graph.sizeVertices();
		this.nEdges = graph.sizeEdges();
		
		this.communities = new int[nVertices];
		for(int i=0; i<nVertices; i++) {
			communities[i] = i;
		}
		
		/* 
		 * matrixE is the matrix of the fractions of edges linking two different communities of the previous 'communities' array.
		 * arrayA is the array of the total fraction of outer edges of the communities of the previous 'communities' array.
		*/
		this.matrixB = computeMatrixB(graph, communities);
		this.arrayK = computeArrayK(matrixB);
		
			
		int[] argmaxModularity = communities.clone();
		double modularity = computeModularity(graph, communities);
		double maxModularity = -2;
		MaxModularityIncreaseParams maxModularityIncreaseParams = new MaxModularityIncreaseParams();
		
		int firstCommunity;
		int secondCommunity;
		double maxModularityIncrease;
		
		for(int i=0; i<(nVertices-1); i++){
			maxModularityIncreaseParams = findMaxModularityIncreaseParams();
			maxModularityIncrease = maxModularityIncreaseParams.maxModularityIncrease;
			firstCommunity = maxModularityIncreaseParams.firstCommunityMaxModularityIncrease;
			secondCommunity = maxModularityIncreaseParams.secondCommunityMaxModularityIncrease;
			
			/* We update the matrix B and the array K using the formulae of the papers. */
			arrayK[firstCommunity] = arrayK[firstCommunity] + arrayK[secondCommunity];
			
			for (int j=0; j<nVertices; j++){
				if (j != firstCommunity) {
					matrixB[firstCommunity][j] = matrixB[firstCommunity][j] + matrixB[secondCommunity][j];
					matrixB[j][firstCommunity] = matrixB[firstCommunity][j];
				} else {
					matrixB[firstCommunity][firstCommunity] = matrixB[firstCommunity][firstCommunity] 
														    + matrixB[secondCommunity][secondCommunity]
															+ 2 * matrixB[firstCommunity][secondCommunity];
				}
			}
			
			/* We merge the community of 'firstCommunity' and 'secondCommunity'.*/
			for (int j=0; j<nVertices; j++) {
				if (communities[j] == secondCommunity) {
					communities[j] = firstCommunity;
				}
			}
			
			modularity += maxModularityIncrease;
			if (modularity > maxModularity) {
				maxModularity = modularity;
				argmaxModularity = communities.clone();
			}
		}
		for(int i=0; i<nVertices; i++) {
			communities[i] = 0;
		}
		endTime = System.nanoTime();
		double duration= (endTime - startTime)/1000000;
		System.out.printf("Greedy algo.: %f milliseconds\n", duration);
		return argmaxModularity;
		
	}

}
