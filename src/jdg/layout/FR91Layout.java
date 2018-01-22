package jdg.layout;

import jdg.graph.AdjacencyListGraph;
import jdg.graph.Node;

import Jcg.geometry.Point_3;
import Jcg.geometry.Vector_3;

/**
 * A class implementing the force directed algorithm by Fruchterman and Reingold (1991)
 * 
 * @author Luca Castelli Aleardi, Ecole Polytechnique
 * @version fev 2017
 */
public class FR91Layout extends Layout {
	// parameters of the algorithm by Fruchterman and Reingold
	public double k; // natural spring length
	public double area; // area of the drawing (width times height)
	public double C; // step
	public double temperature; // initial temperature
	public double minTemperature; // minimal temperature (strictly positive)
	public double coolingConstant; // constant term: the temperature decreases linearly at each iteration
	
	public int iterationCount=0; // count the number of performed iterations
	private int countRepulsive=0; // count the number of computed repulsive forces (to measure time performances)
	
	public int[] indices;
	public int nVertices;
	
	/**
	 * Initialize the parameters of the force-directed layout
	 * 
	 *  @param graph  input graph to draw
	 *  @param width  width of the drawing area
	 *  @param height  height of the drawing area
	 */
	public FR91Layout(AdjacencyListGraph graph, double width, double height) {
		System.out.print("Initializing force-directed method: Fruchterman-Reingold 91...");
		if (graph == null) {
			System.out.println("Input graph not defined");
			System.exit(0);
		}
		
		this.graph = graph;
		this.indices = graph.getIndices();
		this.nVertices = graph.sizeVertices();
		
		// set the parameters of the algorithm FR91
		this.C = 1.;
		this.width = width;
		this.height = height;
		this.area = width * height;
		this.k = C * Math.sqrt(area/nVertices);
		this.temperature = width/2.; // the temperature is a fraction of the width of the drawing area
		this.minTemperature = 0.05;
		this.coolingConstant = 0.99;
		
		System.out.println("done ("+ nVertices +" nodes)");
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
		return (k * k)/distance;
	}

	/**
	 * Perform one iteration of the Force-Directed algorithm.
	 * Positions of vertices are updated according to their mutual attractive and repulsive forces.
	 */	
	public void computeLayout() {
		System.out.print("Performing iteration (FR91): " + this.iterationCount);
		
		/*
		 * For evaluating time performances.
		 */
		long startTime = System.nanoTime();
		long endTime; 

		// First step: for each vertex compute the displacements due to attractive and repulsive forces.
		Vector_3 force;
		double forceNorm;
		Node node;
		Vector_3[] attractiveForces = computeAllAttractiveForces();
		Vector_3[] repulsiveForces = computeAllRepulsiveForces();
		Vector_3[] displacements = new Vector_3[nVertices];
		
		for (int i=0; i<this.nVertices; i++) {
			node = this.graph.getNode(this.indices[i]);
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
		
        // update temperature
        this.cooling(); 
		
		// evaluate time performances
    	endTime=System.nanoTime();
        double duration = (double)(endTime - startTime)/1000000000.;
        System.out.println("iteration " + this.iterationCount + " done (" + duration + " seconds)");
	
		this.iterationCount++; // increase counter (to count the number of performed iterations)
	}
	
	/**
	 * Compute the displacement of vertex 'u', due to repulsive forces (of all nodes)
	 * 
	 * @param u  the vertex to which repulsive forces are applied
	 * @return 'displacement' a 3d vector storing the displacement of vertex 'u'
	 */	
	private Vector_3 computeRepulsiveForce(Node node) {
		Vector_3 force = new Vector_3(0, 0, 0);
		Vector_3 vector;
		Point_3 nodePoint = node.getPoint();
		double forceNorm;
		double scalar;
		for (Point_3 v: this.graph.listOfPoints()) {
			if (!v.equals(nodePoint)) {
				vector = (Vector_3) (nodePoint.minus(v));
				forceNorm = Math.sqrt((double)vector.squaredLength());
				vector = vector.divisionByScalar(forceNorm);
				scalar =  - repulsiveForce(forceNorm);
				force = force.sum(vector.multiplyByScalar(scalar));
			}
		}
		return force;
	}
	
	/**
	 * Compute, for each vertex, the displacement due to repulsive forces (between all nodes)
	 * 
	 * @return a vector v[]: v[i] stores the geometric displacement of the i-th node
	 */	
	private Vector_3[] computeAllRepulsiveForces() {
		Vector_3[] repulsiveForces = new Vector_3[this.nVertices];
		for (int i: this.indices) {
			repulsiveForces[i] = computeRepulsiveForce(this.graph.getNode(i));
		}
		return repulsiveForces;
	}
	
	/**
	 * Compute the displacement of vertex 'u', due to the attractive forces of its neighbors
	 * 
	 * @param u  the vertex to which attractive forces are applied
	 * @return 'disp' a 3d vector storing the displacement of vertex 'u'
	 */	
	private Vector_3 computeAttractiveForce(Node u) {
		Vector_3 force = new Vector_3(0, 0, 0);
		Vector_3 vector;
		double vectorNorm;
		double scalar;
		for (Node v: this.graph.getNeighbors(u)) {
			vector = (Vector_3) (u.getPoint().minus(v.getPoint()));
			vectorNorm = Math.sqrt((double)vector.squaredLength());
			vector = vector.divisionByScalar(vectorNorm);
			scalar = attractiveForce(vectorNorm);
			force = force.sum(vector.multiplyByScalar(scalar));
		}
		return force;
	}
	
	/**
	 * Compute, for each vertex, the displacement due to attractive forces (between neighboring nodes)
	 * 
	 * @return a vector v[]: v[i] stores the geometric displacement of the i-th node
	 */	
	private Vector_3[] computeAllAttractiveForces() {
		Vector_3[] attractiveForces = new Vector_3[this.nVertices];
		for (int i: indices) {
			attractiveForces[i] = computeAttractiveForce(this.graph.getNode(i));
		}
		return attractiveForces;
	}
	
	/**
	 * Cooling system: the temperature decreases linearly at each iteration
	 * 
	 * Remark: the temperature is assumed to remain strictly positive (>=minTemperature)
	 */	
	protected void cooling() {
		this.temperature = Math.max(this.temperature * this.coolingConstant, this.minTemperature);
		//this.temperature=Math.max(this.temperature-coolingConstant, minTemperature); // variant
	}
	
	public String toString() {
		String result = "force-directed algorihm: Fruchterman Reingold\n";
		result = result + "\t area= " + this.width + " x " + this.height+ "\n";
		result = result + "\t k= " + this.k + "\n";
		result = result + "\t C= " + this.C + "\n";
		result = result + "\t initial temperature= " + this.temperature + "\n";
		result = result + "\t minimal temperature= " + this.minTemperature + "\n";
		result = result + "\t cooling constant= " + this.coolingConstant + "\n";
		return result;
	}
	
}
