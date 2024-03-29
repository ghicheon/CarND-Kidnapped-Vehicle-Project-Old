/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"


#define  PI 3.141592 

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

    struct Particle p;
    int i;
    default_random_engine gen;

    //normal (Gaussian) distribution for x,y,theta
    normal_distribution<double> dist_x(x, std[0]);
    normal_distribution<double> dist_y(y, std[1]);
    normal_distribution<double> dist_theta(theta,std[2]);


    num_particles=10;       

//    particles.reserve(100); //make room for 100 particles!
//   weights.reserve(100); //make room for 100 particles!

    p.weight = 1;

    for(i=0; i <  num_particles ; i++)
    {
        p.id = i;
        p.x = dist_x(gen);
        p.y = dist_y(gen);
        p.theta = dist_theta(gen);

        particles.push_back(p);
        weights.push_back(1);
    }

    is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/


    default_random_engine gen;

    for(int i=0; i < particles.size() ; i++ ) 
    {
        double x0 = particles[i].x;
        double y0 = particles[i].y;
        double theta0 = particles[i].theta;

        double x;
        double y;
        double theta;

        if( yaw_rate == 0 )
        {
            x = x0 + velocity*delta_t*cos(theta0);
            y = y0 + velocity*delta_t*sin(theta0);
            theta = theta0;
        }
        else
        {
            x = x0 + velocity/yaw_rate * ( sin(theta0 + yaw_rate * delta_t) - sin(theta0));
            y = y0 + velocity/yaw_rate * ( cos(theta0) - cos(theta0 + yaw_rate * delta_t) );
            theta = theta0 + yaw_rate * delta_t;
        }

        normal_distribution<double> dist_x(x, std_pos[0]);
        normal_distribution<double> dist_y(y, std_pos[1]);
        normal_distribution<double> dist_theta(theta,std_pos[2]);

        particles[i].x = dist_x(gen);
        particles[i].y = dist_y(gen);
        particles[i].theta = dist_theta(gen);
    }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

    /* !!! this operation is implementated in updateWeights() because of efficiency. !!!  */

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

    default_random_engine gen;

//return;


    double finalWeight = 1.0;
    double min=99999;
    double stdxx = std_landmark[0]*std_landmark[0];
    double stdyy = std_landmark[1]*std_landmark[1];
    double stdxy = std_landmark[0]*std_landmark[1];

    for(int i=0; i < particles.size() ; i++ ) 
    {
        double xp = particles[i].x;
        double yp = particles[i].y;
        double theta = particles[i].theta;

        double x;
        double y;

        finalWeight = 1.0;

        for(int ii=0 ; ii < observations.size() ; ii++)
        {
            LandmarkObs obs = observations[ii];
            double xc = obs.x;
            double yc = obs.y;

            //transform observations(car) system  into global map system. 
            double xm =  xp + cos(theta)*xc - sin(theta)*yc;
            double ym =  yp + sin(theta)*xc + cos(theta)*yc;

            int closest_mark = 0 ; 

            min=99999;

            //find closest landmark from this observation!
            for(int k=0; k< map_landmarks.landmark_list.size() ; k++)
            {
                Map::single_landmark_s    mark = map_landmarks.landmark_list[k] ;

                double d = dist( xm , ym , mark.x_f ,mark.y_f );
                if( d < min )
                {
                    min = d;
                    closest_mark = k;
                }
            }

            //set the landmark id to this observation.
            // obs.id = map_landmarks.landmark_list[closest_mark].id_i;

            double x = xm;
            double y = ym;

            double ux= map_landmarks.landmark_list[closest_mark].x_f;
            double uy= map_landmarks.landmark_list[closest_mark].y_f;

            double w  = -(  ( (x-ux)*(x-ux))/(2*stdxx) + ((y-uy)*(y-uy))/(2*stdyy) );

            //calculate weight for this observation.
            double vvv = 1/(2*PI* stdxy) * exp(w);

            finalWeight *=   vvv;

        }
        weights[i] = finalWeight;
        particles[i].weight =  finalWeight;
    }

}



void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    default_random_engine gen;
    discrete_distribution<int> dist(weights.begin(),weights.end());

    vector<Particle> particles_new;


    for(int i=0; i < num_particles; i++)
    {
        particles_new.push_back(particles[dist(gen)]);
    }
    particles = particles_new ;

}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
