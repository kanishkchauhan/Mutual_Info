/*                                  MODIFIED VERSION (USES JUST ONE FUNCTION)


 *        File: firing_VS_coupling.cc
 *      Author: Kanishk Chauhan
 *        Date: June 8, 2019
 * Description: to determine the dependence of the firing rate, the coefficient of variance 
 *              , and count of spikes on the coupling strength and the number of nodes in a 
 *              Hodgkin-Huxley model.
 *               
 */


# include <iostream>
# include <cmath>
# include <cstdlib> //needed for random numbers
# include <random> // for different distributions of random numbers
# include <fstream> // needed to write the threshold value data to a text file

using namespace std;

// following function will update the values of V, m, and h at each tim e step
void update(double &v, double v0[], double &m, double &h, double Iext, double kappa, int N, int k, double dt, double noise);

const int ttime = 50; // total time of simulation in seconds
const int tmax = ttime * 1000; // total time in ms
const double dt = 0.0001; //time step
const double tu = 100; //relaxation time

int main ()
{
    // N is the number of nodes including the central one, the Nth node being the central one
    // kappa is the coupling constant

    //following chunk of codes is needed to increment the value of kappa logarithmically 
    double kappa_min = 0.4;
    double kappa_max = 1000;
    double number_of_points = 20; // number of values of kappa
    double alpha = 1/(number_of_points - 1) * log10(kappa_max/kappa_min);
    // this 'alpha' will be used to obtain the next value of kappa

    ofstream file; //creating the file to store data
    file.open("firing_VS_coupling_data.csv");
    file << "Coupling,F-R central,CV central,Count of Spikes" << endl;
    
    double spike_mat[20000];// to store the spike times for the central node
    double D = 500; // noise intensity
    double dn = sqrt(2*dt*D)/2; 
    double nu, nt, rnd1, rnd2, r, theta, noise;
    double Iext = 20.0; // external current

    // initiallizing random number generator
    default_random_engine generator;
    uniform_real_distribution<double> random(0,1); // now if I write "random (generator)", it gives a uniform random number between 0 and 1

    int n[] = {3,6,10};
    for (int N : n) 
    {
        file << endl;
        file << "For N = " << N << endl;
        double v0[N], m0[N], h0[N]; // arrays for initial conditions
        double variation_amp_v0[N], variation_amp_m0[N], variation_amp_h0[N];
        
        double kappa = kappa_min; // initialization for while loop  
        int counter = 1; //to obtain next value of kappa 

        while (kappa <= kappa_max)
        { 
            double spike_time; 
            int count_cen = 0;// to count the number of spikes in the central node
            counter ++; // this is for incrementing the value of kappa
            
            for (int i =0; i <= N-1; i++)
            {
                // we begin with generating different initial values of V,m,and h for the nodes
                v0[i] = -80; m0[i] = 0.2; h0[i] = 0.5;
                variation_amp_v0[i] = rand()%20-10 ; // small variations in values of v for nodes 
                variation_amp_m0[i] = (rand()%10-5) * 0.01;//small variations in values of m for nodes
                variation_amp_h0[i] = (rand()%10-5) * 0.01;//small variations in values of h for nodes
                v0[i] = v0[i] + variation_amp_v0[i];
                m0[i] = m0[i] + variation_amp_m0[i];
                h0[i] = h0[i] + variation_amp_h0[i];
            } // creates different initial conditions for all nodes by adding small variations
            
            // we start with transition to steady state
            nu = tu/dt + 1; 
            for (double j = 1; j <= nu; j += 1)
            {   
                for (int k = 0; k < N; k++)
                {
                    // generating random numbers for noise term
                    rnd1 = random(generator); // creates one unifirmly distributed random number between 0 and 1
                    rnd2 = random(generator);
                    r = sqrt(-2*log(rnd1));// r and theta are used to create normally distributed numbers from uniformly distributed numbers
                    theta = 2*3.14*rnd2;
                    noise = dn * r * cos(theta);
                    
                    double v = v0[k], m = m0[k], h = h0[k]; // assigns initial values of V,m, and h corresponding to the kth node to new variables to pass through funtions
                    update(v, v0, m, h, Iext, kappa, N, k, dt, noise);
                    v0[k] = v; m0[k] = m; h0[k] = h; // updating the the values
                    
                } // transition to steady ends here for kth node
            } // transition to steady ends here for all N nodes

            nt = tmax/dt;
            double v_new[N]; double m_new[N]; double h_new[N]; //needed for updating the values 
                                                                //of V,m,and h                                                                
            int iflag = 1;
            for (double jj = 1; jj <= nt; jj++) // this updates the values of V,m,and h for each node one by one
            {
                for (int kk = 0; kk < N; kk++)
                {
                    // generating random numbers for noise term
                    rnd1 = random(generator); // creates one unifirmly distributed random number between 0 and 1
                    rnd2 = random(generator);
                    r = sqrt(-2*log(rnd1));// r and theta are used to create normally distributed numbers from uniformly distributed numbers
                    theta = 2*3.14*rnd2;
                    noise = dn * r * cos(theta);
                    
                    double v = v0[kk], m = m0[kk], h = h0[kk]; // assigns initial values of V,m, and h corresponding to the kth node to new variables to pass through funtions
                    update(v, v0, m, h, Iext, kappa, N, kk, dt, noise);
                    v_new[kk] = v; m_new[kk] = m; h_new[kk] = h; // updating the the values
                        
                    // discriminating spikes for the central node
		    
		            if (v0[N-1] > -20 && v_new[N-1] < -20)
                    {
                        iflag = 1;
                    }
                    if (v0[N-1] < 20 && v_new[N-1] > 20 && iflag == 1)
                    {
                        spike_time = (jj - 1)*dt; // this stores spike time for kkth node
                        spike_mat[count_cen] = spike_time;
			            count_cen ++;
                        iflag = 0;
                    }
		    
                    //updating the values of V,m, and h
                    v0[kk] = v_new[kk];
                    m0[kk] = m_new[kk];
                    h0[kk] = h_new[kk];                    
                }    
            }// spike_mat now has all spike time values for the central node 

            // calculations for inter-spike interval and firing rate
            // the following array will store the isi for the central node
            double isi_central[count_cen-1];
            double sum2 = 0, sum1 = 0;
            double meanISI_central, firing_rate_central, CV_central;
            
            for (int val = 0; val < count_cen-1; val++) // this is to find mean inter spike interval
            {
                isi_central[val] = spike_mat[val+1] - spike_mat[val];
                sum1 += isi_central[val];
            }
            meanISI_central = sum1 / (count_cen-1);
            firing_rate_central = 1000/meanISI_central;

            for (int vall = 0; vall < count_cen-1; vall++) // this is for standard deviation in inter spike interval
            {
                double num = pow((isi_central[vall] - meanISI_central),2);
                sum2 += num;
            }
            double sd = sqrt(sum2/(count_cen-1));  // standard deviation
            CV_central = sd/meanISI_central;  // covariance
	        cout << "central firing rate = " << firing_rate_central << " and CV = " << CV_central << endl;
            
            file << kappa << "," << firing_rate_central << "," << CV_central << "," << count_cen << endl; 
            double expo = alpha * (counter - 1);// for updating the value of kappa
            kappa = kappa_min * pow(10,expo);
        }// firing rates, etc for each value of kappa have been calculated and written in a file            
    }
    file.close();
    return 0;
}

void update(double &v, double v0[], double &m, double &h, double Iext, double kappa, int N, int k, double dt, double noise)
{
    double am, bm, ah, bh, fm, fh, fV, I_ion;
    double v1 = v + 20.4; double v2 = v + 25.7; double v3 = v + 114.0;

    am = (1.314 * v1) / (1.0 - exp(-v1/10.3)); // calculating the change in the value of m
    bm = -(0.0608 * v2) / (1.0 - exp(v2/9.16));
    fm = am * (1.0 - m) - bm * m;
    m = m + (dt * fm);

    ah = -(0.068 * v3) / (1.0- exp(v3/11.0)); // calculating the change in the value of h
    bh = 2.52 / (1.0 + exp(-(v+31.8)/13.4));
    fh = ah * (1.0 - h) - bh * h;
    h = h + (dt * fh);

    I_ion = 20.0 * (v + 80.0) + m*m*m*h * 1100.0 * (v - 50.0);
    if (k != N-1)
    {   
        fV = (-I_ion + Iext + kappa * (v0[N-1] - v))/2;
        v = v + (dt * fV) + noise;            
    }
    else
    {
        double sum = 0;
        for (int i = 0; i < N-1; i++)
        {
            sum += v0[i];
        }
        double fV = (-I_ion - kappa * ((N-1)*v - sum))/2;
        v = v + (dt * fV);
    }
}