package jdg.layout;

import java.util.Random;

import jdg.graph.AdjacencyListGraph;
import jdg.graph.Node;

import Jcg.geometry.Point_3;

/**
 * Network visualization: abstract implementation of the force-directed paradigm
 * 
 * @author Luca Castelli Aleardi, Ecole Polytechnique (INF421)
 * @version october 2017
 */
public abstract class Layout {
	public AdjacencyListGraph graph;
	public double width, height; // dimensions of the drawing area
	
	public static int seed = 10;
	/** Random generator */	
	static Random generator = new Random(seed); // initialize random generator
	
	/**
	 * Initialize vertex locations at random (in a square of given size WxH)
	 */	
	public static void setRandomPoints(AdjacencyListGraph graph, double width, double height) {
		Point_3 point;
		double width1 = width/2., height1=height/2.;
		for(Node u: graph.vertices){
			double n1 = width1 - 2 * width1 * generator.nextDouble();
			double n2 = height1 - 2 * height1 * generator.nextDouble();
		    point = new Point_3 (n1, n2, 0.0);
		    u.setPoint(point);
		}		
	}
	
	/**
	 * Perform one iteration of the Force-Directed algorithm.
	 * Positions of vertices are updated according to their
	 * mutual attractive and repulsive forces.
	 */	
	public abstract void computeLayout();
			
}
