/**
 * @file 	problem.c
 * @brief 	Example problem: shearing sheet.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail 	This problem uses shearing sheet boundary
 * conditions. Particle properties resemble those found in
 * Saturn's rings.
 *
 * @section 	LICENSE
 * Copyright (c) 2011 Hanno Rein, Shangfei Liu
 *
 * This file is part of rebound.
 *
 * rebound is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * rebound is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "main.h"
#include "particle.h"
#include "boundaries.h"
#include "output.h"
#include "input.h"
#include "communication_mpi.h"
#include "tree.h"
#include "tools.h"
#include "display.h"
#include "integrator.h"
#include "collisions.h"
#include "collision_resolve.h"


extern double OMEGA;
extern double minimum_collision_velocity;

extern double (*coefficient_of_restitution_for_velocity)(double);
double coefficient_of_restitution_bridges(double v);

extern double opening_angle2;
void problem_damping_force();
void collision_resolve_merger(struct collision c);


void problem_init(int argc, char* argv[]){
	// Setup constants
#ifdef GRAVITY_TREE
	opening_angle2	= .5;
#endif // GRAVITY_TREE
	integrator 			= SEI;
//	integrator			= IAS15;
	collision_resolve 		= collision_resolve_merger;			// Setup our own collision routine.

	boxsize 			= 50;
	OMEGA 				= 1.;	// 1/s
	G 				= 1.;		// N / (1e-5 kg)^2 m^2
	double Q 			= 3.;
	double Tau 			= .01;
	double N 			= 300.;
	double cs 			= Q * M_PI;
	double heating 	  		= Q * Q * M_PI * M_PI / boxsize / boxsize / 8.;
//	softening 			= 20. * sqrt(4. * Tau / ( M_PI * N * cs)) * boxsize;
	softening 			= 2. * sqrt(4. * Tau / (M_PI * N)) * boxsize;
//	softening 			= 0;
	//softening 			= 0.001;			// m
	//dt 				= 1e-4*2.*M_PI/OMEGA;	// s
//	softening 			= 0.488;			// m
	dt 				= 1.e-2 * 2. * M_PI / OMEGA / 6.23;	// s - correction factor brings it into normal units for comparison with estimate 
	tmax				= 1.e2 * 2. * M_PI / OMEGA;
#ifdef OPENGL
	display_rotate_z		= 20;			// Rotate the box by 20 around the z axis, then
	display_rotate_x		= 60;			// rotate the box by 60 degrees around the x axis
#ifdef LIBPNG
	system("mkdir png");
#endif // LIBPNG
#endif // OPENGL
	root_nx = 2; root_ny = 2; root_nz = 1;
	nghostx = 2; nghosty = 2; nghostz = 0; 			// Use two ghost rings
	double surface_density 		= 1.;
//	double particle_density		= 4.108;   //   - this is the density of the individual particles! since M_sun=3.6e14 in sim
	double particle_radius_min 	= softening / 2.;  // softening is a way to avoid the singular 1/r potential problems
//	double particle_radius_min 	= sqrt(Tau / (4. * M_PI * N)) * boxsize;
	double particle_density		= 3. * boxsize * boxsize / (M_PI * N * particle_radius_min * particle_radius_min * particle_radius_min);				   // 1/r -- > 1/sqrt(r^2 + softening^2)
							   // want it to be about 2 - 3 times the radius of the particle - collision of cotton balls
							   // and not lead balls, avoid gravitational scattering - unimportant in reality, as
							   // stuff is mostly very light
//	double ringmassfraction		= 0.2; // fraction of total mass in an artificial ring
//	double ringwidthfraction	= 0.1; //  ring width
//	if (argc>1){			// Try to read boxsize from command line
//		boxsize = atof(argv[1]);
//	}
	init_box();

	// Initial conditions
	// Use Bridges et al coefficient of restitution.
	coefficient_of_restitution_for_velocity = coefficient_of_restitution_bridges;
	//minimum_collision_velocity = particle_radius_min*OMEGA*1e7;  // small fraction of the shear
	minimum_collision_velocity = particle_radius_min*OMEGA*1e2;  // small fraction of the shear
	double total_mass = 4.*surface_density*boxsize*boxsize;
//	double ring_mass  = ringmassfraction*total_mass;
//	double ring_width = ringwidthfraction*boxsize_x;
	double mass       = 0;
//	double heating    = 3.15e-4;   // precursor to speed of sound
//	double heating 	  = Q / boxsize_x * M_PI;
		while(mass < total_mass){
			struct particle pt;
			pt.x 		= tools_uniform(-boxsize_x/2.,boxsize_x/2.);
			pt.y 		= tools_uniform(-boxsize_y/2.,boxsize_y/2.);
			pt.z 		= tools_normal(heating/1.414)*boxsize_x; // normal generates a number about 100 times the variance
			pt.vx 		= tools_normal(heating*2.0)*OMEGA*boxsize_x; //x-direction has 2x
			pt.vy 		= -1.5*pt.x*OMEGA+OMEGA*boxsize_y*tools_normal(heating);
			pt.vz 		= tools_normal(heating/1.414)*boxsize_x*OMEGA;
	//		pt.vz 		= 0;
			pt.ax 		= 0;
			pt.ay 		= 0;
			pt.az 		= 0;
			double radius 	= particle_radius_min;
	#ifndef COLLISIONS_NONE
			pt.r 		= radius;						// m
	#endif
			double particle_mass = particle_density*4./3.*M_PI*radius*radius*radius;
			pt.m		 = particle_mass ;

	/* This block below will make sure that particles are not created overlapping with each other. Presently it is still possible for them to overlap across the boundary. */
			int i = 0;
			int problem = 0;
			if (mass == 0)
			{
				particles_add(pt);
				mass += particle_mass;
				continue;
			}
			if (mass == 1)
			{
				particles_add(pt);
				mass += particle_mass;
				continue;
			}
			else
			{
				while(i < mass / particle_mass)
				{

					double x21  = pt.x - particles[i].x;
					double y21  = pt.y - particles[i].y;
					double z21  = pt.z - particles[i].z;

					if (4.*particle_radius_min*particle_radius_min > x21*x21 + y21*y21 + z21*z21)
					{
						problem += 1;
					}
					i += 1;
				}

				if (problem == 0)
				{
					particles_add(pt);
					mass += particle_mass;
				}
			}
	
}
//	mass = 0;
//	while(mass<ring_mass){
//		struct particle pt;
//		pt.x 		= tools_uniform(-ring_width/2.,ring_width/2.);
//		pt.y 		= tools_uniform(-boxsize_y/2.,boxsize_y/2.);
//		pt.z 		= tools_normal(heating/1.414)*ring_width;
//		pt.vx 		= tools_normal(heating*2.0)*ring_width*OMEGA; //x-direction has 2x
//		pt.vy 		= -1.5*pt.x*OMEGA+ring_width*tools_normal(heating)*OMEGA;
//		pt.vz 		= tools_normal(heating/1.414)*ring_width*OMEGA;
//		pt.ax 		= 0;
//		pt.ay 		= 0;
//		pt.az 		= 0;
//		double radius 	= particle_radius_min;
//		pt.r 		= radius;						// m
//		double particle_mass = particle_density*4./3.*M_PI*radius*radius*radius;
//		pt.m		 = particle_mass ;
//		particles_add(pt);
//		mass += particle_mass;
//	}
//	double cs = heating*OMEGA*boxsize_x; // dispersion velocity, scale height = heating*boxsize - there we go, speed of sound
//	double cs = 2. * sqrt(2. * heating) * boxsize;
	printf("Heating: %f\n", heating);
	printf("Mass: %f\n", particle_radius_min*particle_radius_min*particle_radius_min*4*M_PI*particle_density / 3.);
	printf("Radius: %f\n", particle_radius_min);
	printf("Density: %f\n", particle_density);

	printf("Cs: %f\n", cs);
	printf("Toomre wavelength: %f\n",2.*M_PI*M_PI*surface_density/OMEGA/OMEGA*G);
	printf("Initial Toomre Q: %f\n",cs*OMEGA/M_PI/G/surface_density);
//	printf("Optical depth Tau: %f\n",3.*surface_density*cs/particle_radius_min/particle_density/4.);
	printf("Optical depth Tau: %f\n",3.*surface_density/particle_radius_min/particle_density/4.);
//	printf("Ring wavelength: %f\n",cs*cs/M_PI/G/surfacedensity/1.5);

//	problem_additional_forces = problem_damping_force;  //Set function pointer to add dissipative forces.
}


void collision_resolve_merger(struct collision c){
	struct particle p1 = particles[c.p1];
	struct particle p2 = particles[c.p2];
	double x21  = p1.x  - p2.x;
	double y21  = p1.y  - p2.y;
	double z21  = p1.z  - p2.z;
	double rp   = p1.r+p2.r;

	if (rp*rp < x21*x21 + y21*y21 + z21*z21) return; // not overlapping
	double vx21 = p1.vx - p2.vx;
	double vy21 = p1.vy - p2.vy;
	double vz21 = p1.vz - p2.vz;

	if (vx21*x21 + vy21*y21 + vz21*z21 >0) return; // not approaching

	if (p1.lastcollision>=t || p2.lastcollision>=t) return; // already collided


//	if (vx21*vx21 + vy21*vy21 + vz21*vz21 > 10000) {
//		void collision_resolve_hardsphere(struct collision c);
//		return;
//		} // impact velocity too large
	particles[c.p2].lastcollision = t;
	particles[c.p1].lastcollision = t;


	// Note: We assume only one collision per timestep.
	// Setup new particle (in position of particle p1. Particle p2 will be discarded.
	struct particle cm = tools_get_center_of_mass(p1, p2);
	particles[c.p1].x = cm.x;
	particles[c.p1].y = cm.y;
	particles[c.p1].z = cm.z;
	particles[c.p1].vx = cm.vx;
	particles[c.p1].vy = cm.vy;
	particles[c.p1].vz = cm.vz;
	particles[c.p1].r = p1.r*pow(cm.m/p1.m,1./3.);	// Assume a constant density
	particles[c.p1].m = cm.m;
	// Remove one particle.

	particles[c.p2].x = boxsize *1000.;

	// Make sure we don't drift, so let's go back to the center of momentum
//	tools_move_to_center_of_momentum();
}


double coefficient_of_restitution_bridges(double v){
	// assumes v in units of [m/s]
	//double eps = 0.32*pow(fabs(v)*100.,-0.234);
	//if (eps>1) eps=1;
	//if (eps<0) eps=0;
        double  eps=.5;
	return eps;
}

//void problem_damping_force(){
//        return;   //wyq: return if don't want to use this
//        double c = 1e0;
//        for (int i=0;i<N;i++){
//                particles[i].ax -= c*dt*particles[i].vx;
//                particles[i].ay -= c*dt*(particles[i].vy+1.5*particles[i].x*OMEGA);
//                particles[i].az -= c*dt*particles[i].vz;
//        }
//}


void problem_inloop(){
}

void problem_output(){
/*	if (output_check(0.1*2.*M_PI / OMEGA))*/
/*	{*/
/*		output_binary("restart.bin");*/
/*		printf("\nSaved binary file. Restart simulation with './rebound --restart restart.bin'.\n");*/
/*	} */

#ifdef LIBPNG
	if (output_check(.5*2.*M_PI/OMEGA)){
//		output_png("png/");
	}
#endif //LIBPNG
	if (output_check(1e-6*2.*M_PI/OMEGA)){
		output_timing();
		//output_append_velocity_dispersion("veldisp.txt");
	}
	if (output_check(2.*M_PI/OMEGA)){
//		output_ascii("position.txt");
//		printf("%f", collisions_N)
		output_append_ascii("position.txt");
//		double heating = 0.09424;
//		double numb = tools_normal(heating*2.0)*boxsize;
//		printf("%f", numb);
//		printf("%f", heating);



	}
}


void problem_finish(){
}
