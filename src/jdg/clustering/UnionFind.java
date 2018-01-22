package jdg.clustering;
import java.util.HashMap;

/**
 * Implementation of an UnionFind with an hashtable. 
 * @author stolen :)
 *
 */
public class UnionFind {
    private HashMap<Integer, Integer> parent;
    
    public UnionFind( ){
		this.parent = new HashMap<>();
    }
    
    public Integer find( Integer src ){
    	Integer res = src;
		while (this.parent.get(res) != null){
		    res = this.parent.get(res);
		}
		if (!res.equals(src))
		    this.parent.put(src, res);

		return res;
    } 
    public void union(Integer v0, Integer v1){
    	Integer v0_p = find(v0);
    	Integer v1_p = find(v1);
		if (!v1_p.equals(v0_p))
		    this.parent.put(v1_p, v0_p);
		return;
	
    }
    
}
