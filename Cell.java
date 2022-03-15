import java.util.Scanner;
import java.lang.IllegalArgumentException;
public class Cell {
	private final int x;
	private final int y;
	private Grain g = null;
	private boolean boundary;
	private double gamma = -1.0;
	private double[] d = new double[]{0,0,0,0};
	public double curvature;
	public Grain prevGrain;

	public Cell(int x, int y) {
		this.x = x;
			this.y = y;
	}
	
	public void setD(double[] d) {
		this.d = d;
	}

	public double[] getD() {
		return d;
	}

	public void setD0(){
		for(int i=0; i<d.length; i++)
			this.d[i] = 0;
	}

	public int getX() {
		return this.x;
	}

	public int getY() {
		return this.y;
	}
	
	public Grain getGrain() {
		return this.g;
	}

	public void setGrain(Grain gr) {
		this.prevGrain = this.g;
		if (this.g == null) {
			this.g = gr;
			//  System.out.println("Set to: " + this.g.getGrNum());
			//  gr.n_cells++;
		}
		else if ((this.g != gr) || (this.g == gr)) {
			//  System.out.println("yes");
			//  Scanner scan = new Scanner(System.in);
			//  scan.nextInt();
			this.g.n_cells--;
			this.g = gr;
			//   gr.n_cells++;
		}
		this.g.n_cells++;
	}

	public boolean isBoundary() {
		return this.boundary;
	}

	public void setBoundary(boolean b){
		this.boundary = b;
	}

	public double getGamma(int delta_O) {
		this.gamma = 0.8;
		if(Math.abs((double)delta_O) < 15 && delta_O != 0)
			this.gamma = 0.8 * Math.abs((double)delta_O) * (1 - Math.log(Math.abs((double)delta_O)/15))/15;
		else if(delta_O == 0)
			this.gamma =0;
		if (this.gamma < 0)
			System.out.println("Negative gamma");
		return this.gamma;
	}

	public double getGamma() {
		if (this.gamma < 0)
		  throw new IllegalArgumentException("Gamma not set yet. "+ this.y + " " + this.x + " " + this.gamma);
		return this.gamma;
	}
}
