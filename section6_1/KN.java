package section6_1;
import java.util.ArrayList;
import java.util.Random;

public class KN {
	 private int k;
	 private int n0;
	 private double alpha,delta; 
	 private ArrayList<Double> mu = new ArrayList<Double>();
	 private ArrayList<Double> sigma2 = new ArrayList<Double>();
	 private int bestID = -1;
	 private double totalSampleSize=0d;
	 
	 private Random R= new Random();
		
	 public KN() {
		 n0 = 10;
		 alpha = 0.05;
		 delta = 0.1;
		 mu.add(0.1);mu.add(0d);
		 sigma2.add(1d);sigma2.add(1d);
		 R.setSeed(234234);
	 }
	 
	 public KN(int n0, double delta ,double alpha, ArrayList<Double> mu, ArrayList<Double> sigma2,long seed){
	        this.n0 = n0;
		 	this.alpha = alpha;
		 	this.delta = delta;
	        this.mu.addAll(mu);
	        this.sigma2.addAll(sigma2);
	        R.setSeed(seed);
	 }
	 
	 public int getBestID() {
		 return bestID;
	 }
	 
	 public double getTotalSampleSize() {
		 return totalSampleSize;
	 }
	 
	 public void run() {
		 
	
		 k  = mu.size();
		
		 double eta = 0.5*(Math.pow(2*alpha/(k-1),-2.0/(n0-1))-1);
		 
		 double h2 = 2*eta*(n0-1);
		 ArrayList<Integer> I = new ArrayList<Integer>();
		 for(int i = 0 ; i < k ; i++) {
			 I.add(i);
		 }
		 double[][] Xij =  new double[k][k];
		 double[][] Xij2 = new double[k][k];
		 double[][] Sij2 = new double[k][k];
		 
		 double[] tempX = new double[I.size()];
		 
		 for(int ell = 0; ell < n0; ell++) {
			 for(int i = 0; i < I.size(); i++) {
				 tempX[i]= R.nextGaussian() * Math.sqrt(sigma2.get(I.get(i)))+mu.get(I.get(i));
			 }
			 
			 for(int i = 0; i < I.size(); i++) {
				 Xij[I.get(i)][I.get(i)] = Xij[I.get(i)][I.get(i)] + tempX[i];
				 Xij2[I.get(i)][I.get(i)] = Xij2[I.get(i)][I.get(i)] + tempX[i] * tempX[i];
				 
				 for(int j = i + 1; j < I.size(); j++) {
					 Xij[I.get(i)][I.get(j)] = Xij[I.get(i)][I.get(j)]+tempX[i]-tempX[j];
					 Xij[I.get(j)][I.get(i)] = Xij[I.get(j)][I.get(i)]+tempX[j]-tempX[i];
					 Xij2[I.get(i)][I.get(j)] = Xij2[I.get(i)][I.get(j)] +(tempX[i]-tempX[j])*(tempX[i]-tempX[j]);
					 Xij2[I.get(j)][I.get(i)] = Xij2[I.get(i)][I.get(j)];
				 }
			 }
		 }
		 
		 for(int i = 0; i < I.size(); i++) {
			 for(int j = i + 1; j < I.size(); j++) {
				 Sij2[I.get(i)][I.get(j)] = (Xij2[I.get(i)][I.get(j)] -  Xij[I.get(i)][I.get(j)]* Xij[I.get(i)][I.get(j)]/n0)/(n0-1);
				 Sij2[I.get(j)][I.get(i)] = Sij2[I.get(i)][I.get(j)];
			 }
		 }
		 int t = n0;
		 
		 while(I.size()>1) {
			 for(int i  = 0;  i < I.size(); i++) {
				 for (int j = 0; j < I.size(); j++) {
					 if(i!=j) {
						 double Zij = Xij[I.get(i)][I.get(j)];
						 double gt = h2*Sij2[I.get(i)][I.get(j)]/(2.0*delta)-delta*t/2d;
						 if(Zij > gt) {
							 totalSampleSize = totalSampleSize + t;
							 I.remove(j);
							 if(j < i) {
								 i--;
							 }
							 j--;
						 }
					 }
				 }
			 }
			 double[] tempX2 = new double[I.size()];
			 for(int i = 0; i < I.size(); i++) {
				tempX2[i]= R.nextGaussian() * Math.sqrt(sigma2.get(I.get(i)))+mu.get(I.get(i));
			 }
			 
			 for(int i = 0; i < I.size(); i++) {
				 
				 for(int j = i + 1; j < I.size(); j++) {
					 Xij[I.get(i)][I.get(j)] = Xij[I.get(i)][I.get(j)]+tempX2[i]-tempX2[j];
					 Xij[I.get(j)][I.get(i)] = Xij[I.get(j)][I.get(i)]+tempX2[j]-tempX2[i];
					// if(I.get(i)==0&&I.get(j)==1) {
						// System.out.println(Xij[I.get(i)][I.get(j)]/t+" "+t);
					 //}
				 }
			 }
			 if(I.size()==1) {
				 totalSampleSize = totalSampleSize + t;
			 }
			 t = t+1;
		 }
		 bestID = I.get(0);
	 
	 }
	 
}
