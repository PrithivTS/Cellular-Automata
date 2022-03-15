import java.util.*;

public class Grain {
	private final int grNum;
	public int O;
	public int n_cells;
	public double sF;
	public double P;
	private ArrayList<Cell> boundaryCells = new ArrayList<Cell>();
	private int[] d = new int[4];

	public Grain(int grNum, double P, int O, double sF) {
		this.grNum = grNum;
		this.P = P;
		this.O = O;
		this.sF = sF;
		// boundaryCells.add(c);
	}

	public void addBoundaryCell(Cell c) {
		this.boundaryCells.add(c);
	}
	
	public ArrayList<Cell> getBoundaryCells() {
		return this.boundaryCells;
	}
	
	public int getGrNum(){
	return this.grNum;
	}

	public double getDia(double cellSize) {
		return Math.sqrt(n_cells*Math.pow(cellSize * 1E-6,2)*4/Math.PI);
	}
	
	public void clearBoundaryCells() {
		this.boundaryCells.clear();
	}
}
