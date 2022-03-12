/* Restraint to calculate a Gaussian-score based restraint using the coordinates of an anchor point */

#ifndef IMPCENTRIOLES_SINGLE_POINT_MIN_GAUSSIAN_RESTRAINT_H  // Header guard
#define IMPCENTRIOLES_SINGLE_POINT_MIN_GAUSSIAN_RESTRAINT_H

#include <IMP/centrioles/centrioles_config.h>  //To define the IMP<module>EXPORT and IMP<module>_BEGIN/END_NAMESPACE
#include <IMP/Restraint.h>
#include <IMP/core/XYZ.h>
#include <IMP/algebra/Vector3D.h>

IMPCENTRIOLES_BEGIN_NAMESPACE

class IMPCENTRIOLESEXPORT SinglePointMinGaussianRestraint : public IMP::Restraint {
	double sigma_;  //the sigma for the gaussian
	algebra::Vector3D center_;  //x, y, z coordinates of the anchor point
	IMP::ParticlesTemp plist_;  //all the particles for which the minimum score is calculated
	
	public:
		SinglePointMinGaussianRestraint(IMP::ParticlesTemp plist, algebra::Vector3D center, double sigma);
		
		// unprotected_evaluate calculates the score
		// do_get_inputs returns the particles for which the score was calculated
		virtual double unprotected_evaluate(IMP::DerivativeAccumulator* accum) const IMP_OVERRIDE;
		//IMP_OVERRIDE macro ensures that this overrides (and not overloads) a parent method
		virtual IMP::ModelObjectsTemp do_get_inputs() const IMP_OVERRIDE;
		IMP_OBJECT_METHODS(SinglePointMinGaussianRestraint);  //add the usual IMP object methods
	
	private:
		double get_min_distance() const;
};

IMPCENTRIOLES_END_NAMESPACE

#endif
