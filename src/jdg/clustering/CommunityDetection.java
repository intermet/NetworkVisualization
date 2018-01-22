package jdg.clustering;
import jdg.graph.AdjacencyListGraph;
import jdg.graph.Node;

/**
 * This class provides methods for computing the modularity of a partition and
 * for solving the community detection problem
 * 
 * @author Luca Castelli Aleardi (INF421, 2017)
 */
public abstract class CommunityDetection {

	/**
	 * This method returns a partition of a network of size 'n' into communities. <p>
	 * <p>
	 * Remarks:<p>
	 * -) the nodes of the networks are numbered 0..n-1
	 * -) the graph is partitioned into 'k' communities, which have numbers 0..k-1
	 * 
	 * @param graph  the input network (adjacency list representation)
	 * @return an array of size 'n' storing, for each vertex, the index of its community (a value between 0..,k-1)
	 */
	public abstract int[] computeClusters(AdjacencyListGraph graph);
	
	/**
	 * This method computes the matrix 'B' of the 'graph' whose nodes are regrouped into 'communities'.
	 * The matrix B is defined in the report.
	 * 
	 * @param graph  the input network (adjacency list representation).
	 * @param communities an array of size 'n' storing, for each vertex, the index of its community (a value between 0..,k-1).
	 * @return matrixB the matrix B.
	 */
	public int[][] computeMatrixB(AdjacencyListGraph graph, int[] communities){
		int nCommunities = -1;
		for (int i: communities)
			nCommunities = Math.max(nCommunities, i) + 1;
		int matrixB[][] = new int[nCommunities][nCommunities];
		int nodeCommunity;
		int neighbourCommunity;
		for (Node node: graph.vertices) {
			nodeCommunity = communities[node.index];
			for (Node neighbour: node.neighbors) {
				neighbourCommunity = communities[neighbour.index];
				matrixB[nodeCommunity][neighbourCommunity] += 1;
			}
		}
		return matrixB;
	}
	
	
	/**
	 *  This method computes the array 'K' of the graph using its 'matrixB'. This array is defined in the report. 
	 * @param matrixB the matrix B of the graph.
	 * @return arrayK	the array 'K'.
	 */
	public int[] computeArrayK(int[][] matrixB) {
		int nCommunities = matrixB.length;
		int arrayK [] = new int[nCommunities];
		for(int i=0; i<nCommunities; i++) {
			for (int j=0; j<nCommunities; j++) {
				arrayK[i] += matrixB[i][j];
			}
		}
		return arrayK;
	}
	
	/**
	 * This method computes the modularity of a partition of the input graph whose nodes are regrouped into 'k' communities.
	 * We use the formula of the report, using the matrix A and the array K.
	 * Remarks:
	 * -) the nodes of the networks are numbered 0..n-1
	 * -) the graph is partitioned into 'k' communities, which have numbers 0..k-1
	 * 
	 * @param graph  the input network (adjacency list representation)
	 * @param communities  an array of size 'n' storing, for each vertex, the index of its community (a value between 0..,k-1)
	 * @return modularity  the modularity of the partition
	 */
	public double computeModularity(AdjacencyListGraph graph, int[] communities) {
		int nEdges = graph.sizeEdges();
		int[][] matrixB = computeMatrixB(graph, communities);
		int[] arrayK = computeArrayK(matrixB);
		int nCommunities = matrixB.length;
		double modularity = 0;
		for (int i=0; i<nCommunities; i++) {
			modularity += matrixB[i][i] - (arrayK[i] * arrayK[i])/(2. * nEdges);
		}
		return modularity/(2 * nEdges);
	}
	
}
