package section6_1;
import java.util.ArrayList;
import java.util.Random;

import org.apache.commons.math3.distribution.TDistribution;

public class FHN1 {
	 private int k;
	 private int n0;
	 private double alpha; 
	 private ArrayList<Double> mu = new ArrayList<Double>();
	 private ArrayList<Double> sigma2 = new ArrayList<Double>();
	 private int bestID = -1;
	 private double totalSampleSize=0d;
	
	 private Random R= new Random();
	
	 public FHN1() {
		 n0 = 10;
		 alpha = 0.05;
		 mu.add(0.1);mu.add(0d);
		 sigma2.add(1d);sigma2.add(1d);
		 R.setSeed(234234);
	 }
	 public FHN1(int n0, double alpha, ArrayList<Double> mu, ArrayList<Double> sigma2,long seed){
	        this.n0 = n0;
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
		 double[][] Sij2 = new double[k][k];
		 double c;
		 if(n0 == 10) {
			 c = 11.4;//12.6;//14.6;//7.358;
		 }else if(n0==15){
			 c = 12.6;//5.870;
		 }else if(n0==20){
			 c = 13.4;
		 }else if(n0==25){
			 c = 5.040;
		 }else if(n0==30){
			 c = 14.6;//4.671;
		 }else if(n0==35) {
			 c = 4.514;
		 }else {
			 c = 4.415;
		 }
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
						 double Zij = Xij[I.get(i)][I.get(j)]/Sij2[I.get(i)][I.get(j)];
						 double tau = t/Sij2[I.get(i)][I.get(j)];
						 if(Zij > ht(tau,c)) {
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
				 tempX2[i]= R.nextGaussian() * Math.sqrt(sigma2.get(i))+mu.get(i);
			 }
			 
			 for(int i = 0; i < I.size(); i++) {
				 Xij[I.get(i)][I.get(i)] = Xij[I.get(i)][I.get(i)] + tempX2[i];
				 Xij2[I.get(i)][I.get(i)] = Xij2[I.get(i)][I.get(i)] + tempX2[i] * tempX2[i];
				 
				 for(int j = i + 1; j < I.size(); j++) {
					 Xij[I.get(i)][I.get(j)] = Xij[I.get(i)][I.get(j)]+tempX2[i]-tempX2[j];
					 Xij[I.get(j)][I.get(i)] = Xij[I.get(j)][I.get(i)]+tempX2[j]-tempX2[i];
				 }
			 }
			 if(I.size()==1) {
				 totalSampleSize = totalSampleSize + t;
			 }
			 t = t+1;
		 }
		 bestID = I.get(0);
	 }
	 public double ht(double t,double c) {
		 
		 double  ht =Math.sqrt((c+Math.log(t+1))*(t+1));
		 
		 return ht;
	 }
}
