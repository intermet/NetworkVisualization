package jdg.clustering;

import java.util.HashMap;

import jdg.graph.AdjacencyListGraph;
import jdg.graph.Node;

public class LouvainAlgorithm extends CommunityDetection {
	
	public AdjacencyListGraph graph;
	public int nVertices;
	public int nEdges;
	public int[] communities;
	public int[][] matrixSigma;
	public int[][] weights;
	public int[] arrayL;
	public int[] arrayK;
	public UnionFind union;
	

	public LouvainAlgorithm() {
	}
	
	/**
	 * This method computes the matrix \Sigma defined in the report.
	 * @return matrixSigma
	 */
	public int[][] computeMatrixSigma(){
		int[][] matrixSigma = new int[this.nVertices][this.nVertices];
		for (int i=0; i<this.nVertices; i++) {
			for (int j=0; j<this.nVertices; j++) {
				matrixSigma[i][j] = this.weights[i][j];
			}
		}
		return matrixSigma;
	}
	

	 /**
	  * This method computes the array L defined in the report.
	  * @return arrayL
	  */
	public int[] computeArrayL() {
		int[] arrayL = new int[this.nVertices];
		for (int i=0; i<this.nVertices; i++) {
			for (int j=0; j<this.nVertices; j++) {
				arrayL[i] += this.weights[i][j];
			}
		}
		return arrayL;
	}
	
	
	/**
	 * This methods implements the first step of Louvain's algorithm.
	 */
	public void computeLocalMaxima() {
		boolean noChange = false;
		while (!noChange) {
			noChange = increaseModularity();
		}
		return;
	}
	
	
	/**
	 * This methods tries to improve modularity by changing the community of a node to the community of one of its neighbors.
	 * @return
	 */
	public boolean increaseModularity() {
		boolean noChange = true;
		double maxModularityIncrease = -2;
		Node argMaxModularityIncrease = new Node(0);
		double modularityIncrease;
		for (Node node: this.graph.vertices) {
			maxModularityIncrease = -2;
			for (Node neighbour: node.neighbors) {
				if (this.communities[node.index] != this.communities[neighbour.index]) {
					modularityIncrease = computeModularityIncrease(node, neighbour);
					if (modularityIncrease > maxModularityIncrease) {
						maxModularityIncrease = modularityIncrease;
						argMaxModularityIncrease = neighbour;
					}
				}
			}
			if (maxModularityIncrease > 0) {
				changeCommunity(node, argMaxModularityIncrease);
				noChange = false;
			}
		}
		return noChange;
	}
	
	/**
	 * This methods computes the increase of modularity by moving 'node' to the community of 'neighbour'.
	 * It uses the formula given in the report.
	 * @param node
	 * @param neighbour
	 * @return modularityIncrease
	 */
	public double computeModularityIncrease(Node node, Node neighbour) {
		double modularityIncrease = 0;
		int nodeIndex = node.index;

		int nodeCommunity = this.communities[nodeIndex];
		int neighbourCommunity = this.communities[neighbour.index];
		int kNode = this.arrayK[nodeIndex];
		modularityIncrease = this.matrixSigma[nodeIndex][neighbourCommunity] - this.matrixSigma[nodeIndex][nodeCommunity]
						   + kNode*(this.arrayL[nodeCommunity] - this.arrayL[neighbourCommunity])/(2.*this.nEdges)
						   +(this.weights[nodeIndex][nodeIndex] - (kNode*kNode)/(2. * this.nEdges));
		modularityIncrease = modularityIncrease/(this.nEdges);
		return modularityIncrease;
	}
	
	
	/**
	 * This methods moves 'node' to the community of 'neighbour', and updates the matrix \Sigma and the array L
	 * according to formulae given in the report.
	 * @param node
	 * @param neighbour
	 */
	public void changeCommunity(Node node, Node neighbour) {
		int nodeIndex = node.index;
		int nodeCommunity = this.communities[nodeIndex];
		int neighbourIndex = neighbour.index;
		int neighbourCommunity = this.communities[neighbourIndex];
		int nodeNeighbourIndex;
		this.arrayL[nodeCommunity] -= this.arrayK[nodeIndex];
		this.arrayL[neighbourCommunity] += this.arrayK[nodeIndex];
		for (Node nodeNeighbour: this.graph.vertices) {
			nodeNeighbourIndex = nodeNeighbour.index;
			this.matrixSigma[nodeNeighbourIndex][nodeCommunity] -= this.weights[node.index][nodeNeighbour.index];
			this.matrixSigma[nodeNeighbourIndex][neighbourCommunity] += this.weights[node.index][nodeNeighbour.index];
		}
		this.communities[nodeIndex] = this.communities[neighbourIndex];
		return;
	}
	
	/**
	 * This method achieves the second step of Louvain's algorithm by aggregating the graph according to communities.
	 * It also computes the weights of the edges of the graph.
	 * @return
	 */
	public AdjacencyListGraph aggregateGraph() {
		AdjacencyListGraph aggregatedGraph = new AdjacencyListGraph();
		HashMap<Integer, Node> communityMapNode = new HashMap<>();
		Node aggregatedGraphNode;
		Node nodeRepresentant;
		int nodeRepresentantIndex;
		int nodeCommunity;
		int nodeIndex;
		int index = 0;
		
		/* We create one node per community.*/
		for (Node node: this.graph.vertices) {
			nodeIndex = node.index;
			nodeCommunity = this.communities[nodeIndex];
			if (communityMapNode.containsKey(nodeCommunity)) {
				nodeRepresentant = communityMapNode.get(nodeCommunity);
				this.union.union(nodeRepresentant.numericLabel, node.numericLabel);
			} else {
				aggregatedGraphNode = new Node(index);
				aggregatedGraphNode.numericLabel = node.numericLabel;
				aggregatedGraph.addNode(aggregatedGraphNode);
				communityMapNode.put(nodeCommunity, aggregatedGraphNode);
				index += 1;
			}
		}
		
		int[][] aggregatedGraphWeights = new int[index][index];
		
		Node neighbourRepresentant;
		int neighbourRepresentantIndex;
		int neighbourIndex;
		int neighbourCommunity;
		
		/* We computes the new weights */
		for (Node node: this.graph.vertices) {
			nodeIndex = node.index;
			nodeCommunity = this.communities[nodeIndex];
			nodeRepresentant = communityMapNode.get(nodeCommunity);
			nodeRepresentantIndex = nodeRepresentant.index;
			aggregatedGraphWeights[nodeRepresentantIndex][nodeRepresentantIndex] += this.weights[nodeIndex][nodeIndex];
			for (Node neighbour: node.neighbors) {
				neighbourIndex = neighbour.index;
				neighbourCommunity = this.communities[neighbourIndex];
				neighbourRepresentant = communityMapNode.get(neighbourCommunity);
				neighbourRepresentantIndex = neighbourRepresentant.index;
				aggregatedGraph.addEdge(nodeRepresentant, neighbourRepresentant);
				aggregatedGraphWeights[nodeRepresentantIndex][neighbourRepresentantIndex] += this.weights[nodeIndex][neighbourIndex];
			}
		}

		this.weights = aggregatedGraphWeights;
		return aggregatedGraph;
	}
		
	/**
	 * This method achieves the Louvain's algorithm.
	 * @return finalCommunities
	 */
	public int[] computeClusters(AdjacencyListGraph graph) {
		long startTime = System.nanoTime();
		long endTime;
		this.graph = graph;
		this.nVertices = graph.sizeVertices();
		this.nEdges = graph.sizeEdges();
		this.weights = new int[this.nVertices][this.nVertices];
		for (Node node: graph.vertices) {
			node.numericLabel = node.index;
			for (Node neighbour: node.neighbors) {
				this.weights[node.index][neighbour.index] = 1;
			}
		}
		this.matrixSigma = computeMatrixSigma();
		this.arrayL = computeArrayL();
		this.arrayK = computeArrayL();
		this.communities = new int[this.nVertices];
		for (int i=0; i<this.nVertices; i++) {
			this.communities[i] = i;
		}
		
		/* That union lets us save the communities computed during the passes.*/
		this.union = new UnionFind();
		

		boolean noChange = false;
		AdjacencyListGraph currentGraph;
		while(!noChange) {
			computeLocalMaxima();
			currentGraph = aggregateGraph();
			if (currentGraph.sizeVertices() == this.nVertices) {
				noChange = true;
			} else {
				this.graph = currentGraph;
				this.nVertices = this.graph.sizeVertices();
				this.matrixSigma = computeMatrixSigma();
				this.arrayL = computeArrayL();
				this.arrayK = computeArrayL();
				this.communities = new int[this.nVertices];
				for (int i=0; i<this.nVertices; i++) {
					this.communities[i] = i;
				}
			}
		}
		endTime = System.nanoTime();
		double duration= (endTime - startTime)/1000000;
		System.out.printf("Louvain's algo.: %f milliseconds\n", duration);

		int[] finalCommunities = new int[graph.sizeVertices()];
		for (int i=0; i<graph.sizeVertices(); i++) {
			finalCommunities[i] = this.union.find(i);
		}
		return finalCommunities;
	}

}
