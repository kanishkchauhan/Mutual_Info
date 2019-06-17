/*
 *        File: firingVScoupling.cc
 *      Author: Kanishk Chauhan
 *        Date: June 8, 2019
 * Description: to determine the dependence of the firing rate, the coefficient of variance 
 *              , and count of spikes on the intensty of noise in a  Hodgkin-Huxley model.
 *               
 */

# include <iostream>
# include <cmath>
# include <cstdlib> //needed for random numbers
# include <random> // for different distributions of random numbers
# include <fstream> // needed to write the threshold value data to a text file

using namespace std;

// following function will update the values of V, m, and h at each tim e step
void update(double &v, double &m, double &h, double I_ext, double dt, double noise);

const int ttime = 10; // total time of simulation in seconds
const int tmax = ttime * 1000; // total time in ms
const double dt = 0.0001; //time step
const double tu = 100; //relaxation time

int main ()
{
    //following chunk of codes is needed to increment the value of D (noise intensity) logarithmically 
    double Dmin = 10;
    double Dmax = 10000;
    double number_of_points = 20; // number of values of D
    double alpha = 1/(number_of_points - 1) * log10(Dmax/Dmin);
    // this 'alpha' will be used to obtain the next value of D

    ofstream file; //creating the file to store data
    file.open("firing_VS_noise.csv");
    file << "noise intensity (D),Firing Rate,CV,Count of Spikes" << endl;
    
    double spike_mat[20000];// to store the spike times for the central node
    double nu, nt, noise, v, v0, m, m0, h, h0, dn, spike_time, rnd1, rnd2, r, theta;

    // initiallizing random number generator
    default_random_engine generator;
    uniform_real_distribution<double> random(0,1); // now if I write "random (generator)", it gives a unform random number between 0 and 1

    int Iext[] = {20,25,30}; // these values of external current will be used 
    for (int I : Iext) 
    {
        file << endl;
        file << "For Iext = " << I << endl;
        
        int counter = 1; // needed to obtain next value of D 
        double D = Dmin; // initialization for while loop
        while (D < Dmax)
        { 
            v0 = -80; m0 = 0.2; h0 = 0.5;
            dn = sqrt(2*dt*D)/2; 
            counter ++; // this is for incrementing the value of D
            
            // we start with transition to steady state
            nu = tu/dt + 1; 
            for (double j = 1; j <= nu; j += 1)
            {   
                rnd1 = random(generator); // creates one unifirmly distributed random number between 0 and 1
                rnd2 = random(generator);
                r = sqrt(-2*log(rnd1));// r and theta are used to create normally distributed numbers from uniformly distributed numbers
                theta = 2*3.14*rnd2;
                noise = dn * r * cos(theta);
                update(v0, m0, h0, I, dt, noise);
            } // transition to steady ends here

            nt = tmax/dt;
            double v_new; double m_new; double h_new; //needed for updating the values 
                                                                //of V,m,and h    

            int count = 0;// to count the number of spikes in the node                                                            
            int iflag = 1;
            for (double jj = 1; jj <= nt; jj++) // this updates the values of V,m,and h for each node one by one
            {
                rnd1 = random(generator); // creates one unifirmly distributed random number between 0 and 1
                rnd2 = random(generator);
                r = sqrt(-2*log(rnd1));// r and theta are used to create normally distributed numbers from uniformly distributed numbers
                theta = 2*3.14*rnd2;
                noise = dn * r * cos(theta);
                v = v0; m = m0; h = h0;
                update(v, m, h, I, dt, noise);
                v_new = v; m_new = m; h_new = h; // updating the the values    

                // discriminating spikes for the central node
		        if (v0 > -20 && v_new < -20)
                {
                    iflag = 1;
                }
                if (v0 < 20 && v_new > 20 && iflag == 1)
                {
                    spike_time = (jj - 1)*dt; // this stores spike time
                    spike_mat[count] = spike_time; // count here stores the number of spikes
			        count ++;
                    iflag = 0;
                }
                v0 = v_new; m0 = m_new; h0 = h_new; // updating the values for next loop                   
            }// spike_mat now has all spike time values for the central node 

            // calculations for inter-spike interval and firing rate
            // the following array will store the isi for the central node
            double isi[count-1];
            double sum2 = 0, sum1 = 0;
            double meanISI, firing_rate, CV;
            
            for (int val = 0; val < count-1; val++) // this is to find mean inter spike interval
            {
                isi[val] = spike_mat[val+1] - spike_mat[val];
                sum1 += isi[val];
            }
            meanISI = sum1 / (count-1);
            firing_rate = 1000/meanISI;

            for (int vall = 0; vall < count-1; vall++) // this is for standard deviation in inter spike interval
            {
                double num = pow((isi[vall] - meanISI),2);
                sum2 += num;
            }
            double sd = sqrt(sum2/(count-1));  // standard deviation
            CV = sd/meanISI;  // coeff. of variaton
	        cout << "firing rate = " << firing_rate << " and CV = " << CV << endl;
            
            file << D << "," << firing_rate << "," << CV << "," << count << endl; 
            double expo = alpha * (counter - 1);// for updating the value of D
            D = Dmin * pow(10,expo);
        }// firing rates, etc for each value of kappa have been calculated and written in a file            
    }
    file.close();
    return 0;
}

void update(double &v, double &m, double &h, double I_ext, double dt, double noise)
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
    fV = (-I_ion + I_ext)/2;
    v = v + (dt * fV) + noise; 
}