package jdg.layout;

import Jcg.geometry.Point_3;
import Jcg.geometry.Vector_3;
import jdg.graph.AdjacencyListGraph;
import jdg.graph.Node;

/**
 * A class implementing the Fruchterman and Reingold method with fast approximatino of repulsive forces (using octrees)
 * 
 * @author Luca Castelli Aleardi, Ecole Polytechnique
 * @version fev 2017
 */
public class FastFR91Layout extends Layout {
	// parameters of the algorithm by Fruchterman and Reingold
	public double k; // natural spring length
	public double area; // area of the drawing (width times height)
	public double C; // step
	public double temperature; // initial temperature
	public double minTemperature; // minimal temperature (strictly positive)
	public double coolingConstant; // constant term: the temperature decreases linearly at each iteration
	public boolean useCooling; // say whether performing simulated annealing
	
	public int iterationCount = 0; // count the number of performed iterations
	private int countRepulsive = 0; // count the number of computed repulsive forces (to measure time performances)
	
	public int nVertices;
	public int[] indices;
	
	public Octree octree;
	
	/**
	 * Initialize the parameters of the force-directed layout
	 * 
	 *  @param g  input graph to draw
	 *  @param w  width of the drawing area
	 *  @param h  height of the drawing area
	 *  @param C  step length
	 */
	public FastFR91Layout(AdjacencyListGraph graph, double width, double height) {
		System.out.print("Initializing force-directed method: fast Fruchterman-Reingold 91...");
		if(graph == null) {
			System.out.println("Input graph not defined");
			System.exit(0);
		}
		this.graph = graph;
		this.nVertices = graph.sizeVertices();
		this.indices = graph.getIndices();
		
		// set the parameters of the algorithm FR91
		this.C = 1.;
		this.width = width;
		this.height = height;
		this.area = width * height;
		this.k = C * Math.sqrt(area/this.nVertices);
		this.temperature = width/5.; // the temperature is a fraction of the width of the drawing area
		this.minTemperature = 0.05;
		this.coolingConstant = 0.98;
		
		this.octree = new Octree();
		System.out.println("Hello");
		
		System.out.println("done (" + nVertices + " nodes)");
		//System.out.println("k="+k+" - temperature="+temperature);
		System.out.println(this.toString());
	}
	
	/**
	 * Compute the (intensity of the) attractive force between two nodes at a given distance
	 * 
	 * @param distance  distance between two nodes
	 */	
	public double attractiveForce(double distance) {
		return (distance * distance)/k;
	}
	
	/**
	 * Compute the (intensity of the) repulsive force between two nodes at a given distance
	 * 
	 * @param distance  distance between two nodes
	 */	
	public double repulsiveForce(double distance) {
		countRepulsive++;
		return (this.k * this.k)/distance;
	}

	/**
	 * Compute the displacement of vertex 'u', due to the attractive forces of its neighbors
	 * 
	 * @param u  the vertex to which attractive forces are applied
	 * @return 'disp' a 3d vector storing the displacement of vertex 'u'
	 */	
	private Vector_3 computeAttractiveForce(Node node) {
		Vector_3 attractiveForce = new Vector_3(0, 0, 0);
		Vector_3 vector;
		double vectorNorm;
		double scalar;
		for (Node neighbour: this.graph.getNeighbors(node)) {
			vector = (Vector_3) (node.getPoint().minus(neighbour.getPoint()));
			vectorNorm = Math.sqrt((double)vector.squaredLength());
			if (vectorNorm > 0) {
				vector = vector.divisionByScalar(vectorNorm);
				scalar = attractiveForce(vectorNorm);
				attractiveForce = attractiveForce.sum(vector.multiplyByScalar(scalar));
			}
			
		}
		return attractiveForce;
	}
	
	/**
	 * Compute, for each vertex, the displacement due to attractive forces (between neighboring nodes)
	 * 
	 * @return a vector v[]: v[i] stores the geometric displacement of the i-th node
	 */	
	private Vector_3[] computeAllAttractiveForces() {
		Vector_3[] attractiveForces = new Vector_3[this.nVertices];
		for (int i: this.indices) {
			attractiveForces[i] = computeAttractiveForce(this.graph.getNode(i));
		}
		return attractiveForces;
	}	
	
	/**
	 * Compute the displacement of vertex 'u', due to repulsive forces (of all nodes)
	 * 
	 * @param u  the vertex to which repulsive forces are applied
	 * @return 'displacement' a 3d vector storing the displacement of vertex 'u'
	 */	
	private Vector_3 computeRepulsiveForce(Node node, Octree octree) {
		if(octree == null)
			System.out.println("OCTREE NULL");
		Vector_3 repulsiveForce = new Vector_3(0, 0, 0);
		Point_3 nodePoint = node.getPoint();
		
		Vector_3 diff;
		double norm;
		double s;
		Point_3 barycenter = octree.barycenter;
		norm = Math.sqrt((double) (nodePoint.minus(barycenter)).squaredLength());
		// if the node is far from the nodes in octree we use the approximation
		if (octree.squareWidth <= octree.theta * norm) {
			diff = (Vector_3) (node.getPoint().minus(barycenter));
			norm = Math.sqrt((double)diff.squaredLength());
			diff = diff.divisionByScalar(norm);
			s =  - repulsiveForce(norm);
			repulsiveForce = repulsiveForce.sum(diff.multiplyByScalar(s));
			return repulsiveForce;
		} else {
			if (octree.lowleftOctree.nPoints > 0)
				repulsiveForce.sum(computeRepulsiveForce(node, octree.lowleftOctree));
			if (octree.lowrightOctree.nPoints > 0)
				repulsiveForce.sum(computeRepulsiveForce(node, octree.lowrightOctree));
			if (octree.upleftOctree.nPoints > 0)
				repulsiveForce.sum(computeRepulsiveForce(node, octree.upleftOctree));
			if (octree.uprigthOctree.nPoints > 0)
				repulsiveForce.sum(computeRepulsiveForce(node, octree.uprigthOctree));
			return repulsiveForce;
		}
	}
	
	
	/**
	 * Compute, for each vertex, the displacement due to repulsive forces (between all nodes)
	 * 
	 * @return a vector v[]: v[i] stores the geometric displacement of the i-th node
	 */	
	private Vector_3[] computeAllRepulsiveForces() {
		Vector_3[] repulsiveForces = new Vector_3[this.nVertices];
		for (int i: this.indices) {
			repulsiveForces[i] = computeRepulsiveForce(this.graph.getNode(i), this.octree);
		}
		return repulsiveForces;
	}
	

	/**
	 * Perform one iteration of the Force-Directed algorithm.
	 * Positions of vertices are updated according to their mutual attractive and repulsive forces.
	 */	
	public void computeLayout() {
		System.out.print("Performing iteration (fast FR91): " + this.iterationCount);
		
		/*
		 * for evaluating time performances
		 */
		long startTime = System.nanoTime();
		long endTime;
		
		
		this.octree.init(this.octree, graph.listOfPoints(), 0);
		
		Vector_3 force;
		double forceNorm;
		Node node;
		Vector_3[] attractiveForces = computeAllAttractiveForces();
		Vector_3[] repulsiveForces = computeAllRepulsiveForces();
		Vector_3[] displacements = new Vector_3[this.nVertices];
		
		for (int i=0; i<this.nVertices; i++) {
			node = this.graph.getNode(indices[i]);
			force = attractiveForces[i].sum(repulsiveForces[i]);
			forceNorm = Math.sqrt((double) force.squaredLength());
			force = force.divisionByScalar(forceNorm);
			displacements[i] = force.multiplyByScalar(Math.min(this.temperature, forceNorm));
		}
		
		// second step: compute the total displacements and move all nodes to their new locations
        for (int i=0; i<this.nVertices; i++) {
        	node = this.graph.getNode(this.indices[i]);
        	node.setPoint(node.getPoint().sum(displacements[i]));
        }
		
		
		// evaluate time performances
    	endTime = System.nanoTime();
        double duration= (endTime - startTime)/1000000;
        System.out.println("iteration " + this.iterationCount + " done (" + duration + " milliseconds)");

		this.cooling(); // update temperature
	
		this.iterationCount++; // increase counter (to count the number of performed iterations)
	}
	
	/**
	 * Cooling system: the temperature decreases linearly at each iteration
	 * 
	 * Remark: the temperature is assumed to remain strictly positive (>=minTemperature)
	 */	
	protected void cooling() {
		this.temperature = Math.max(this.temperature * coolingConstant, minTemperature);
		//this.temperature=Math.max(this.temperature-coolingConstant, minTemperature); // variant
	}
	
	public String toString() {
		String result = "fast implementation of the force-directed algorihm: Fruchterman Reingold\n";
		result = result +"\t area= " + this.width + " x " + this.height + "\n";
		result = result + "\t k= " + this.k + "\n";
		result = result + "\t C= " + this.C + "\n";
		result = result + "\t initial temperature= " + this.temperature + "\n";
		result = result + "\t minimal temperature= " + this.minTemperature + "\n";
		result = result + "\t cooling constant= " + this.coolingConstant + "\n";
		
		return result;
	}

}
