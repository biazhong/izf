package section6_1;
import java.util.ArrayList;
import java.util.Random;

import org.apache.commons.math3.distribution.TDistribution;


public class procedure2 {
	 private int k; 
	 private double alpha; 
	 private ArrayList<Double> mu = new ArrayList<Double>();
	 private ArrayList<Double> sigma2 = new ArrayList<Double>();
	 private int bestID = -1;
	 private double totalSampleSize=0d;
	 private Random R= new Random();
	 
	 
	 
	 public procedure2() {
		 alpha = 0.05;
		 mu.add(0.1);mu.add(0d);
		 sigma2.add(1d);sigma2.add(1d);
	 }
	 public procedure2(double alpha, ArrayList<Double> mu, ArrayList<Double> sigma2,long seed){
	        this.alpha = alpha;
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
		 ArrayList<Integer> I = new ArrayList<Integer>();
		 for(int i = 0 ; i < k ; i++) {
			 I.add(i);
		 }
		 double[][] Xij =  new double[k][k];
		 double[][] Xij2 = new double[k][k];
		 int r  = 1;
		 while(I.size()>1) {
			 int t = (int)Math.pow(2, r);
			 
			 int n =  t/2;
			 
			 if(r == 1) {
				 n = 2;
			 }
			 
			 
			 double[] tempX = new double[I.size()];
			 
			 for(int ell = 0; ell < n; ell++) {
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
			 
			 for(int i  = 0;  i < I.size(); i++) {
				 for (int j = 0; j < I.size(); j++) {
					 if(i!=j) {
						 double Sij2 = (Xij2[I.get(i)][I.get(j)] -  Xij[I.get(i)][I.get(j)]* Xij[I.get(i)][I.get(j)]/t)/(t-1);
						 double Zij = Xij[I.get(i)][I.get(j)]/Math.sqrt(Sij2);
						 //System.out.println(t+" "+Zij+" "+ht(t)+" "+Math.sqrt(Sij2)+" "+Xij[I.get(i)][I.get(j)]);
						 if(Zij > ht(t)) {
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
			 if(I.size()==1) {
				 totalSampleSize = totalSampleSize + t;
				// System.out.println(t);
			 }
			 
			 r++;
		 }
		 bestID = I.get(0);
		 
		 
	 }
	 
	 
	 public double ht(int t) {
		 
		 double p =  1-alpha/((k-1)*Math.pow(Math.log(2*t)/Math.log(2),2));
		 
		 TDistribution td = new TDistribution(t-1);
		 
		 double  ht = Math.sqrt(t)*td.inverseCumulativeProbability(p);
		 
		 return ht;
	 }
}
