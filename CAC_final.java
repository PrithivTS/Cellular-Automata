import java.awt.*;
import java.awt.image.BufferedImage;
import java.util.*;
import java.io.*;
import java.text.*;
import javax.imageio.ImageIO;
import javax.swing.JFrame;
public class CAC_final
{
	public static void main(String[] args)
	{
		long start=System.currentTimeMillis();
		DecimalFormat df=new DecimalFormat("#.###");
		Scanner scan=new Scanner(System.in);

		Color C[]=new Color[29];
		C[0]=new Color(0,0,0);				//black
		C[1]=new Color(255,0,0);			//red
		C[2]=new Color(0,255,0);			//lime
		C[3]=new Color(0,0,255);			//blue
		C[4]=new Color(255,255,0);			//yellow
		C[5]=new Color(0,255,255);			//aqua
		C[6]=new Color(255,0,255);			//magenta
		C[7]=new Color(192,192,192);			//silver
		C[8]=new Color(128,0,0);			//maroon
		C[9]=new Color(128,128,0);			//olive
		C[10]=new Color(0,128,0);			//green
		C[11]=new Color(128,0,128);			//purple
		C[12]=new Color(0,128,128);			//teal
		C[13]=new Color(0,0,128);			//navy
		C[14]=new Color(255,192,203);			//pink
		C[15]=new Color(255,20,147);			//deep pink
		C[16]=new Color(255,215,0);			//golden
		C[17]=new Color(220,20,60);			//crimson
		C[18]=new Color(240,128,128);			//light coral
		C[19]=new Color(255,165,0);			//orange red
		C[20]=new Color(70,130,180);			//steel blue
		C[21]=new Color(0,191,255);			//deep sky blue
		C[22]=new Color(135,206,235);			//sky blue
		C[23]=new Color(75,61,139);			//indigo
		C[24]=new Color(148,0,211);			//deep violet
		C[25]=new Color(238,130,238);			//violet
		C[26]=new Color(210,105,30);			//chocolate
		C[27]=new Color(230,230,250);			//lavender
		C[28]=new Color(255,255,255); //white
		Random rand=new Random();
		int width =100;
		int height =100;
		double cellSize = 1;				// in micro meter
		int in_dia =13;
		int i, j;
		Cell N[][] = new Cell[height][width];
		/* Initializing all cells */
		for (j = 0; j < width; j++) {
		  for (i = 0; i < height; i++) {
			N[i][j] = new Cell(j, i);
		  }
		}
   
		ArrayList<Grain> grains = new ArrayList<Grain>();
		
		inimicro o = new inimicro(width, height, in_dia, N, rand, grains,cellSize);
		o.generate_microstructure();

		//initializing equations: P_bar - Average dislocation density,
		int percentComplete = 1;
		double burger_vector=2.5E-10, t=0, PAvg = 0, PAvg1 = 0, PMax =0, sr, T,T0, dt,   kb=1.38064852E-23;
		double Totalt=0;
		System.out.print("Enter the strain rate: ");
		sr = scan.nextDouble();//INPUT STRAIN RATE
		System.out.print("Enter the temperature in degree C: " );
		T0 = scan.nextDouble() + 273;//INPUT TEMPERATURE
		double shear_mod=46.4E9;//change shear modulous/ or use relation variable with temperature. 
		double alpha;
		alpha=0.5;
		double curStrain = 0;
		double reqStrain=1.7;//INPUT APPLIED STRAIN
		double beta=1;
		double delStrain = reqStrain / 1000;
		double k1;
		double M;
		double nDot=0;		
		double nucl;
		double gamma=1;
		double curTime=0;
		double reqTime=1.7/sr+0;
		int nNuclei;
		double stressval=0;
		int n_grains = o.n_grains;    
		int n_grains_before_DRX = o.n_grains;	//no. of inital grains
		double avgDia = 0;	//average diameter of grain
		double drxavgDia = 0;	//avergae diameter of grain		
		int[][] directions = new int[][]{{-1,0}, {0,1}, {1,0}, {0, -1}};
		ArrayList<Cell> allBoundaryCells = new ArrayList<Cell>();
		int rNum;
		System.out.println("Number of grains in inital microstructure " + o.n_grains);
		BufferedImage image = new BufferedImage(width,height,BufferedImage.TYPE_INT_RGB);
		int imageNum = 0;
		int count_cells_DRX = 0;
		try {
			PrintStream fs = new PrintStream(new File("flowstress.xls"));
			PrintStream gs = new PrintStream(new File("grainsize.xls"));
			fs.println("Strain\tFlow stress (MPa)");
			gs.println("Strain\tAverage Diameter (micron)\tDRX grain size (micron)\tDRX fraction");
			fs.println((curStrain) + "\t" + (alpha*1E-6*shear_mod*burger_vector*Math.sqrt(PAvg)));
			double DRXF=(double)count_cells_DRX/(width*height);
			gs.println((curStrain) + "\t" + avgDia*1E6 + "\t" + drxavgDia*1E6 + "\t" + DRXF);
			int ab=0;				
			while (curTime < reqTime) {
				int counter = 0; 
				avgDia = 0;
				drxavgDia = 0;
				PMax = 0;
				PAvg = 0;
				PAvg1=0;
				//variable strain rate starts// TO COMPARE WITH FEM. IN FEM, MIDDLE PORTION WILL HAVE HIGH STRAIN RATE FOLLOWING THE RELATION
				double sre = (-1.7*Math.pow(curStrain,2) + 3.3*curStrain + 0.9)*sr;// 1273-0.001
				sre=sr;//IF THE STRAIN RATE IS FIXED. AFTER FEM IMPLEMENTATION LOCAL VARIATION IN STRAIN RATE IS EXPECTED
				// variable strain rate ends
				//
				//adiabatic heating starts
				double ak=0;
				if (sr>=1 & sr < 10){
					ak=211.38615386578-0.13384615386*T0;
				}
				if (sr>=10){
					ak =242.29077-0.14308*T0;
				}			
				
				T =ak*curStrain+T0;
				if (curStrain>reqStrain){
				T =ak*reqStrain+T0;
				}
				//adiabatic heating ends
				//
				//time step begin
				dt =  delStrain/sre;
				//time step ends
				//
				//Zener hollomon begin
				double Qact=482319.715203914;//ACTIVATION ENERGY FOR DEFORMATION
				double Z=sre*Math.exp(Qact/(8.314*T));
				double lnZ=Math.log(Z);
				long Az=9722017481736390L;
				double lnaz=Math. log(Z/Az);
				//Zener hollomon ends
				//
				//k2 begin
				double k2=Math. exp(8.32-0.11*lnZ);
				//k2 ends
				//
				// k1, peak, saturation and yiled stress, steady state, critical  definition
				double Sp;
				Sp=1.05*100*Math.log(Math.pow(Z/Az,1/3.1)+Math.pow((Math.pow(Z/Az,2/3.1)+1),0.5));
				double Yieldstress=Sp/4;
				double Pys=Math.pow((Yieldstress*1E6*2/(burger_vector*shear_mod)),2);
				k1=k2*(Sp)*1E6/(alpha*shear_mod*burger_vector);	
				double theta0=alpha*shear_mod*burger_vector*k1/2;
				double mval=alpha*shear_mod/(k1*in_dia);
				double nfact=0.1;
				double dislocation_energy = alpha*shear_mod * Math.pow(burger_vector,2);
				double Pcr=Math.pow((0.84*k1/k2),2);			// Critical dislocation density
				double Pm=Math.pow((1*k1/k2),2);// MAXIMUM DISLOCATION DENSITY
				//STEADY STATE DEFENITION
				double sfact;
				if (lnZ<40){
					sfact=-0.0527*lnZ + 2.738;
				}
				else{
					sfact=0.0345*lnZ- 0.6481;
				}	
				double Pss=Math.pow((sfact*k1/k2),2);	//steady state stress, sfact is the peak steady state ratio
				// k1, peak, saturation and yiled stress, steady state, critical  definition
				//
				//mobility starts
				double Dob=(6.1E-4)*5E-10;// change lattice diffusivity. 
				double Qfact=239;
				if (T>1323)
				{
				Qfact=239-0.35*(T-1323);	
				} 
				if (T<=1323)
				{ 
				Qfact=239;
					if (T<=1250){						
					Qfact=263;
					}
				}
				if (Qfact<100){
					Qfact=100;
				}
				//moblity ends
				//boundary pressure modification parameter bigins// neglect possible
				if (sre>0.1){
				double c1=0.5*Math.exp(40000/(8.314*T))/Math.exp(40000/(8.314*1273));
				if (c1>0.5){
				c1=0.5;
				}
				beta=Math.pow((sre/0.1),c1);	
				}
				// boundary modification parameter neglected.
				if (beta>2){
					beta=2;
				}
				//boundary pressure modification parameter end
				M=beta*(Dob*burger_vector/(kb*(T)))*Math. exp(-Qfact*1E3/(8.314*(T)));
				//boundary pressure modification parameter end
				//nucleation rate starts
				double fc=0.64;
				if(T>1373){
				fc=0.64+0.016*(T-1373)*0.2;
				}
				if (fc >0.8){
					fc=0.8;
				}
				double nt=-1+4.8*Math.exp(-102*(Math.pow(((T/1223)-1.09),2)));// follow a normal distribution
				double ntl=-1+4.8*Math.exp(-102*(Math.pow(((1223/1223)-1.09),2)));
				if (T>1423){
						 nt=-1+4.8*Math.exp(-102*(Math.pow(((1423/1223)-1.09),2)));
					}
				if (nt<=ntl){
						 nt=ntl;
				}	
				nDot=2E9*nt*Math. pow(sr,fc)/Math. pow(0.001,fc);
				//nucleation rate ends	
				//
				/*
				System.out.println("peakstress "+Sp);
				System.out.println("Ln (A/Z): "+lnaz);
				System.out.println("k1: "+k1 ); 
				System.out.println("k2: "+k2 );
				System.out.println("PS: "+Sp );
				Pss=Math.pow((sfact*k1/k2),2);
				System.out.println("sr " + sr);
				System.out.println("T " + T);
				System.out.println("Qmob " + Qfact);	
				System.out.println("Nucleation rate: " + nDot);	
				System.out.println("Saturation stress: " + Pss);
				System.out.println("Mobility: " + M);
				*/
				//

				
				int curGrains = 0;
				int countDrxGrains = 0;
				allBoundaryCells.clear();
				for (Grain g: grains){
					g.P*=g.sF;
					double A;
					A=g.P;
					double B=0;
					if (curStrain<=reqStrain){
					A =  g.P + ((k1*Math.sqrt(g.P) - k2 * g.P+(1/(in_dia*burger_vector))) * delStrain) ;
					}
				 if(g.getGrNum() > n_grains_before_DRX)
					{	
						if(g.n_cells > 1) 		//Minimum DRX grain - no. of cells
						{
							drxavgDia += g.getDia(cellSize);
							countDrxGrains++;
						}
						if(g.n_cells <= 2){
							g.P=1;
						}
						if(g.n_cells > 2) {
							
							double m=g.getDia(cellSize);
							if (m<0.001E-6){
							m=13E-6;
							}
							if (curStrain<=reqStrain){
							B =  g.P + ((sfact*k1 * Math.sqrt(g.P) - 1*k2 * g.P+(1/(m*burger_vector))) * delStrain) ;
							g.P=Math.min(A, B);
							}
							else {
							if (B>=Pss){
								B=Pss;
							}
							}
						}
					}
					if (B==0){
						g.P=A;
							if(g.P<Pys){
							g.P=Pys;							
							}
					}
					if (B==0 & curStrain>reqStrain){
						g.P=A;
							if(g.P<Pss){
							g.P=Pss;							
							}
					}
					if(g.P > Pcr){
						allBoundaryCells.addAll(g.getBoundaryCells());
					}
					PMax = Math.max(PMax, g.P);
					if(g.n_cells != 0){
						PAvg += g.P * g.n_cells ;
						//PAvg += g.P;
						if(g.n_cells > 0)  //min. no of cells can be changed - normal grain size calc
						{
						avgDia += g.getDia(cellSize);
						curGrains++;
						}
					}					
				}			
		
				avgDia /= curGrains;
				if(countDrxGrains != 0)
					drxavgDia /= countDrxGrains;
				else
					drxavgDia = 0;
				PAvg /= (width*height);
				//System.out.println(allBoundaryCells.size());
				int Np = allBoundaryCells.size();
				//System.out.println(Np);
				if (Np > 0) {
					nNuclei = (int)Math.round((nDot * Np * Math.pow(cellSize * 1E-6,2) * dt));
					//   System.out.println("Nucleating... " + nNuclei);
					 //System.out.println(nNuclei);
					while (nNuclei > 0) {
						Cell c = allBoundaryCells.get(rand.nextInt(Np));
						if ((c.getGrain().P/PMax) > rand.nextDouble()) {    //TODO : Take care of boundary cell updation
						// System.out.println("Nucleation");
						Grain nucleatedGrain = new Grain(++n_grains, 1, 1 + rand.nextInt(180),(0.5+rand.nextDouble())*0+1);
						grains.add(nucleatedGrain); // TODO: Take care that this addition to grains doesn't mess up iteration
						c.setGrain(nucleatedGrain);
						nNuclei--;
						//nucleatedGrain.setCenter(c);		//TODO: Find curvature of new boundary cells and define boundaryCells
						double d[] = new double[]{0,0,0,0};
						c.setD(d);			//consider about marking this as a boundary cell
						//	System.out.println("Nucleation " + c.getX() + " " + c.getY());scan.nextInt();
						}
					}
				}
				// o.markBoundary();

				for(Grain g: grains) {
					ArrayList<BoundaryMigration> coords = new ArrayList<BoundaryMigration>();
					double nd=g.getDia(cellSize);
				if (nd<0.000000000001E-6){
					nd=0.000000000001E-6;
				}
					for(Cell c:g.getBoundaryCells()) {
						double d[] = new double[4];
						int y = c.getY();
						int x = c.getX();
						if (y-1 >= 0)
							d[0] = (M * (dislocation_energy * (gamma)*(N[y-1][x].getGrain().P-g.P) - c.curvature * c.getGamma(g.O-N[y-1][x].getGrain().O)) * dt /(cellSize * 1E-6));
						else
							d[0] = 0;
						if (x-1 >= 0)
							d[1] = (M * (dislocation_energy * (gamma)* (N[y][x-1].getGrain().P-g.P) - c.curvature * c.getGamma(g.O-N[y][x-1].getGrain().O)) * dt /(cellSize * 1E-6));
						else
							d[1] = 0;
						if (x+1 <= width)
							d[2] = (M * (dislocation_energy * (gamma)*(N[y][x+1].getGrain().P-g.P) - c.curvature * c.getGamma(g.O-N[y][x+1].getGrain().O)) * dt /(cellSize * 1E-6));
						else
							d[2] = 0;
						if (y+1 <= height)
							d[3] = (M * (dislocation_energy * (gamma)*(N[y+1][x].getGrain().P-g.P) - c.curvature * c.getGamma(g.O-N[y+1][x].getGrain().O)) * dt /(cellSize * 1E-6));
						else
							d[3] = 0;
						//	if(Np > 0 && x == 32 && y == 762)
						//	System.out.println(Arrays.toString(d));
						double d1[] = c.getD();
						int d2[] = new int[4];
						for (int k=0; k<=3; k++) {
							d[k] += d1[k];
							if(d[k] < 0)
								d[k] = 0;
							else if(d[k] >= 1 ) {
								d2[k] = 1;
								d[k] -= 1;
							}
							else
								d2[k] = 0;
						}
						c.setD(d);
				
						for (i=-d2[0]; i<=d2[3]; i++) {
							for (j=-d2[1]; j<=d2[2]; j++) {
								if((i==0 && j==0) || j+x < 0 || j+x >= width || i+y < 0 || i+y >= height)
									continue;
								if(i==0 || j==0){
									if(0.8 >= rand.nextDouble())
										coords.add(new BoundaryMigration(j+x, i+y, g.getGrNum(),c));
								}
								else{
									if(0.25 >= rand.nextDouble())
										coords.add(new BoundaryMigration(j+x, i+y, g.getGrNum(),c));
								}
							}
						}
					}

					for(BoundaryMigration bm : coords) {
						N[bm.y][bm.x].setGrain(grains.get(bm.grNum-1));
						N[bm.y][bm.x].setD(bm.c.getD());
						bm.c.setD0();
					}
				}
				// System.out.println();
				// System.out.println("Current Strain: "+ df.format(curStrain + delStrain) + "\t Grain size: " + avgDia);
				count_cells_DRX = 0;
				for(Grain g:grains)
					if(g.getGrNum() > n_grains_before_DRX)
						count_cells_DRX += g.n_cells;
				o.markBoundary();
				curStrain += delStrain;
				curTime += dt;
				fs.println((curStrain) + "\t" + (alpha*1E-6*shear_mod*burger_vector*Math.sqrt(PAvg)));
				gs.println((curStrain) + "\t" + avgDia*1E6 + "\t" + drxavgDia*1E6 + "\t" + (double)count_cells_DRX/(width*height));
				// TODO: For all the images check setRGB
	
				if((int)(curTime/reqTime*100) == percentComplete * 10)
				{
					System.out.println(percentComplete*10 + "% complete");
					percentComplete++;
				}
			/*	if(percentComplete == 5)
				{
					System.out.println(curGrains);
					for(Grain g:grains)
						System.out.println(g.getGrNum() + "\t" + g.n_cells);
					scan.nextInt();
				}	*/
				if(imageNum%10 == 0)
				{
					try
					{
						for(i=0; i<height; i++)
						{
							for(j=0; j<width; j++)
							{
								if(N[i][j].isBoundary())
									image.setRGB(j,i,C[0].getRGB());
								else if(N[i][j].getGrain().getGrNum() <= o.n_grains)
									image.setRGB(j,i,C[28].getRGB());
								else {
									int f = N[i][j].getGrain().getGrNum() % 27;
									image.setRGB(j,i,C[f].getRGB());
								}

							}

						}
					String s = imageNum + ".jpg";
					File ouptut = new File(s);
					ImageIO.write(image, "jpg", ouptut);
					}catch(Exception e){System.out.println(e);}
				}
				imageNum++;
				Totalt=Totalt+dt;
			}

			count_cells_DRX = 0;
			for(i=0; i<height; i++)
			{
				for(j=0; j<width; j++)
				{
					if(N[i][j].getGrain().getGrNum() > n_grains_before_DRX)
						count_cells_DRX++;
				}
			}
			System.out.println("Percent DRX is:" + df.format((double)count_cells_DRX/(width*height)));
			//System.out.println("Zener:" + Math.log(Z));
			//System.out.println("k1:" + k1);
			//System.out.println("k2:" + k2);
			//System.out.println("Stress:" + Sp);
			//System.out.println("M:" + M);
			long end=System.currentTimeMillis();
			System.out.println("Computation time: "+ df.format((double)((end-start)/(1000.*60)))+" minutes");
			System.out.println("DRXGS " + drxavgDia*1E6);
			fs.close();
			gs.close();
		}catch(IOException e){System.out.println("File error");}
	}
}
 class inimicro {
	private final int width;
	private final int height;
	private final int in_dia;
	private Random rand;
	private final int n_pixels;
	public int n_grains;
	private Cell N[][];
	private ArrayList<Grain> grains;
	private ArrayList<TriplePoint> coords = new ArrayList<TriplePoint>();
	BufferedImage image;

	public inimicro(int width, int height, int in_dia, Cell[][] N, Random rand, ArrayList<Grain> grains, double cellSize) {
		this.width = width;
		this.height = height;
		this.in_dia = in_dia;
		this.n_pixels = width*height;
		this.n_grains = (int) Math.round(((this.n_pixels * Math.pow(cellSize,2))/(Math.PI*(in_dia*in_dia)/4)));
		this.N = N;
		this.rand = rand;
		this.grains = grains;
		this.image=new BufferedImage(this.width,this.height,BufferedImage.TYPE_INT_RGB);
	}

	public void generate_microstructure() {
		this.nucleate();
		this.grow();
		this.markBoundary();
		try
		{
			for(int i=0; i<this.height; i++)
			{
				for(int j=0; j<this.width; j++)
				{
				Color newcolor=new Color(255,255,255);
				image.setRGB(j,i,newcolor.getRGB());
				if(N[i][j].isBoundary()) {
				  newcolor=new Color(0,0,0);
				  image.setRGB(j,i,newcolor.getRGB());
				}
				}
			}
			File ouptut = new File("init_micro.jpg");
			ImageIO.write(image, "jpg", ouptut);
		}catch(Exception e){System.out.println(e);}
	}

	public void nucleate() {
		int x, y;
		//To store the coordinates of randomly generated nuclei
		int coord_x[]= new int[this.n_grains];
		int coord_y[]= new int[this.n_grains];
		// Coordinates of the first nucleus
		coord_x[0] = this.rand.nextInt(width);
		coord_y[0] = this.rand.nextInt(height);
		// minallowabledist is for keeping sufficent distance between grains to satisfy grain size requirement.
		double minallowabledist = 0.5 * this.in_dia;
		double min_dist = Integer.MAX_VALUE;
		double dist;
		// 1E6*(1+this.rand.nextInt(10) is to keep dislocation density between 1E6 and 1E7.
		// 1 + this.rand.nextInt(180) to get the orientation between 1 and 18
		Grain gr = new Grain(1, 1E6*(1+this.rand.nextInt(10)), 1 + this.rand.nextInt(180),(0.5+rand.nextDouble())*0+1);
		this.N[coord_y[0]][coord_x[0]].setGrain(gr);
		this.grains.add(gr);
		int counter = 1;
		// Nucleating rest of the grains while honoring the grain dia constraint.
		while (counter < this.n_grains) {
			x = this.rand.nextInt(this.width);
			y = this.rand.nextInt(this.height);
			for (int j = 0; j < counter; j++) {
				dist = Math.sqrt((x-coord_x[j])*(x-coord_x[j])+(y-coord_y[j])*(y-coord_y[j]));
				min_dist = Math.min(dist, min_dist);
		}
		if (min_dist >= minallowabledist) {
			// System.out.println("YO");
			coord_x[counter] = x;
			coord_y[counter] = y;
			gr = new Grain(counter+1, 1E6*(1+this.rand.nextInt(10)), 1 + this.rand.nextInt(180),(0.5+rand.nextDouble())*0+1);
			this.N[y][x].setGrain(gr);
			this.grains.add(gr);
			counter++;
		}
		min_dist = Integer.MAX_VALUE;
		}
		System.out.println("Successful nucleation");
	}
	//	Scanner scan = new Scanner(System.in);
	//int fun=0;
	public void grow() {
		// coords contains the cells and their to be marked before next iteration.
		this.coords.clear();
		int i, j, grNum, cx, cy, counter=0;
		double prob;
		// Using Moore type neighborhood
		int[][] directions = new int[][]{{-1,-1}, {-1,0}, {-1,1},  {0,1}, {1,1},  {1,0},  {1,-1},  {0, -1}};

		for (i = 0; i < this.height; i++) {
			for (j = 0; j < this.width; j++) {
				if (this.N[i][j].getGrain() != null){
					grNum = this.N[i][j].getGrain().getGrNum();
					counter = 0;
					for (int[] direction : directions) {
						cx = j + direction[0];
						cy = i + direction[1];
						// prob=0.8 for the 4 adjacent cells and 0.25 for the 4 diagonal cells.
						prob = new double[]{0.8, 0.25}[(++counter)%2];
						double rNum = this.rand.nextDouble();
						// Checking for 1. Validity of coordinates, 2. The cell we are taking doesn't belong to any other Grain
						// and 3. the probability condition is satisfied.
						if (cx < 0 || cx >= this.width || cy < 0 || cy >= this.height || this.N[cy][cx].getGrain() != null || rNum >= prob)
						continue;
						// Marking the cells
						this.coords.add(new TriplePoint(cx, cy, grNum));
					}
				}
			}
		}

		for(TriplePoint tp : this.coords) {
		this.N[tp.y][tp.x].setGrain(this.grains.get(tp.grNum-1));
		}
		// coords's size will be 0, when all the cells are marked. Until then, we grow.
		if (this.coords.size() != 0)
			this.grow();
		else
			this.correctNeighbors();
	}

	public void markBoundary() {
		// directions for curvature calculation by kink method
		int[][] directions = new int[][]{{-2,-2}, {-2,-1}, {-2,0},  {-2,1}, {-2,2},  {-1,2},  {0,2},  {1, 2}, {2, 2}, {2, 1}, {2, 0}, {2, -1}, {2, -2}, {1, -2}, {0, -2}, {-1, -2}, {-1,-1}, {-1,0}, {-1,1},  {0,1}, {1,1},  {1,0},  {1,-1},  {0, -1}};
		for (Grain g: grains) {
			for (Cell c: g.getBoundaryCells())
			c.setBoundary(false);
			g.clearBoundaryCells();
		}

		for (int j = 0; j < width-1; j++) {
			for (int i = 0; i < height-1; i++) {
				//Checking for the right neighbor.
				if (N[i][j].getGrain() != N[i+1][j].getGrain()) {
					N[i][j].setBoundary(true);
					N[i+1][j].setBoundary(true);
				}
				// Checkiong for the bottom neighbor
				if (N[i][j].getGrain() != N[i][j+1].getGrain()) {
					N[i][j].setBoundary(true);
					N[i][j+1].setBoundary(true);
				}
				if (N[i][j].isBoundary())
					N[i][j].getGrain().addBoundaryCell(N[i][j]);
			}
		}

		for (Grain g: grains) {
			for (Cell c: g.getBoundaryCells()){
				// Calculating curvature for boundary cells
				// kink stores the number of cells of the same grain in 25 neighbor system.
				if (c.isBoundary()){
					int kink=1;
					int j = c.getX();
					int i = c.getY();
					for(int[] direction:directions){
						int cx = j + direction[0];
						int cy = i + direction[1];
						if (cx < 0 || cx >= this.width || cy < 0 || cy >= this.height || N[cy][cx].getGrain()!=c.getGrain())
						  continue;
						kink++;
					}
					c.curvature = (1.25 * (15 - kink)/(1E-6 * 25));
				}
			}
		}
	}

	public void correctNeighbors() {
		// Because of the way we are marking the cells, situation may arise where
		// a cell is surrounded at all sides by cells of a single other grain.
		// Then we mark our cell to be part of that grain
		for (int j = 1; j < this.width-1; j++) {
			for (int i = 1; i < this.height-1; i++) {
				if (N[i][j+1].getGrain()==N[i][j-1].getGrain() && N[i][j-1].getGrain()==N[i+1][j].getGrain() && N[i+1][j].getGrain()==N[i-1][j].getGrain()) {
					Grain g = N[i][j+1].getGrain();
					if (N[i][j].getGrain() != g)
					N[i][j].setGrain(g);
				}
			}
		}
	}
}

class TriplePoint {
	public int x;
	public int y;
	public int grNum;
	public TriplePoint(int x, int y, int grNum) {
		this.x = x;
		this.y = y;
		this.grNum = grNum;
	}
}

class BoundaryMigration{
	public int x;
	public int y;
	public int grNum;
	public Cell c;
	public BoundaryMigration(int x, int y, int grNum, Cell c){
		this.x = x;
		this.y = y;
		this.grNum = grNum;
		this.c = c;
	}
}

