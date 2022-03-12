#include <IMP/core/XYZ.h>
#include <IMP/algebra/distance.h>
#include <IMP/algebra/VectorD.h>
#include <IMP/algebra/Vector3D.h>
#include <IMP/centrioles/SinglePointMinGaussianRestraint.h>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <limits>

IMPCENTRIOLES_BEGIN_NAMESPACE

SinglePointMinGaussianRestraint::SinglePointMinGaussianRestraint(IMP::ParticlesTemp plist, algebra::Vector3D center, double sigma) :
		Restraint(plist[0]->get_model(), "SinglePointMinGaussianRestraint %1%"),
		sigma_(sigma),
		center_(center),
		plist_(plist) {}

double SinglePointMinGaussianRestraint::get_min_distance() const {
	double min_distance =  std::numeric_limits<double>::infinity();                         //returns the positive infinity value of the given floating-point type
	for (long unsigned int i = 0; i < plist_.size(); i++){  //find the minimum distance of the particle from the center
		//take the coordinates using the decorator XYZ
		double current_distance = algebra::get_squared_distance(core::XYZ(plist_[i]).get_coordinates(), center_);   //get_coordinates returns a vector with the X, Y, Z coordinates
		if (current_distance <= min_distance) {                                                                     //get_distance returns the distance between two vectors (points)
			min_distance = current_distance;                                                                    //currentMin = currVal; if currVal is closer to the center than currentMin
		}
	}
	return min_distance;
}

double SinglePointMinGaussianRestraint::unprotected_evaluate(IMP::DerivativeAccumulator* accum) const {
	double min_dist = get_min_distance();
	double score = 1. / (sqrt(2. * IMP::PI) * sigma_) * exp(-(min_dist) / (2. * sigma_ * sigma_));      
	
	if (accum){};
	return -log(score);
}

IMP::ModelObjectsTemp SinglePointMinGaussianRestraint::do_get_inputs() const {
	return plist_;
}


IMPCENTRIOLES_END_NAMESPACE
