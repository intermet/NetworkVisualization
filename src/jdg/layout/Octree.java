package jdg.layout;
import java.util.ArrayList;
import java.util.Arrays;

import Jcg.geometry.GeometricOperations_3;
import Jcg.geometry.PointCloud_3;
import Jcg.geometry.Point_;
import Jcg.geometry.Point_3;
import jdg.graph.AdjacencyListGraph;

/**
 * This class implements the structure of octree with trees.
 * @author ziyed
 *
 */
public class Octree {
	
	ArrayList<Point_3> points; // the points described by the octree
	Point_3 barycenter; // the barycenter of 'points'
	int nPoints; //  the number of points
	double squareWidth; // the width of the square where are the points
	double theta = 1.2; // the theta parameter accorded to the paper [3]
	Octree lowleftOctree;
	Octree lowrightOctree;
	Octree upleftOctree;
	Octree uprigthOctree;
	int maxTreeLevel = 11; // the maximum height of the octree
		
	public Octree(){
		this.points = null;
		this.lowleftOctree = null;
		this.lowrightOctree = null;
		this.upleftOctree = null;
		this.uprigthOctree = null;
	}
	
	
	
	public void init(Octree octree, ArrayList<Point_3> points, int level) {
		octree.nPoints = points.size();
		/* We compute the barycenter of the points */
		octree.barycenter = new Point_3();
		Point_3[] pointsArray = new Point_3[points.size()];
		int i = 0;
		for (Point_3 p:points){
			pointsArray[i] = p;
			i += 1;
		}
		octree.barycenter.barycenter(pointsArray);
		
		/**
		 * We stop the propagation of the tree if there is at most one point or if the octree is already to deep.
		 */
		if (points.size() <= 1 || level >= this.maxTreeLevel) {
			octree.points = points;
			return;
		}
		
		// We find the points in the different quarters.
		ArrayList<Point_3> lowleft = new ArrayList<>();
		ArrayList<Point_3> lowright = new ArrayList<>();
		ArrayList<Point_3> upleft = new ArrayList<>();
		ArrayList<Point_3> upright = new ArrayList<>();
		
		// We computes the coordinates of the extremal points
		double minX, minY, maxX, maxY;
		minX = points.get(0).x;
		maxX = minX;
		minY = points.get(0).y;
		maxY = minY;
		for(Point_3 p: points) {
			minX = Math.min(minX, p.x);
			maxX = Math.max(maxX, p.x);
			minY = Math.min(minY, p.y);
			maxY = Math.max(maxY, p.y);
		}
		
		this.squareWidth = Math.max(maxX - minX, maxY - minY); // the width of the square where are the points
		double mid_x = minX + (this.squareWidth/2); // the middle of the square in the direction x
		double mid_y = minY + (this.squareWidth/2); // the middle of the square in the direction y
		
		for (Point_3 p: points) {
			if (p.x <= mid_x && p.y<= mid_y)
				lowleft.add(p);
			else if (p.x <= mid_x)
				upleft.add(p);
			else if (p.y <= mid_y)
				lowright.add(p);
			else
				upright.add(p);
		}
		
		octree.lowleftOctree = new Octree();
		octree.lowrightOctree = new Octree();
		octree.upleftOctree = new Octree();
		octree.uprigthOctree = new Octree();
		// We build the rest of the octree recursively
		init(octree.lowleftOctree, lowleft, level + 1);
		init(octree.lowrightOctree, lowright, level + 1);
		init(octree.upleftOctree, upleft, level + 1);
		init(octree.uprigthOctree, upright, level + 1);
	}
	
	
}
