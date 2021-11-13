#ifndef MONITORINGSTAGEDFRROTREEGENERATORLOGGER_H
#define MONITORINGSTAGEDFRROTREEGENERATORLOGGER_H
#include"../core/GeneratorData.h"
#include"../structures/tree/AbstractCostEstimator.h"
#include"../structures/tree/VolumetricCostEstimator.h"
#include"../structures/tree/SproutingVolumetricCostEstimator.h"
#include"../structures/tree/AdimSproutingVolumetricCostEstimator.h"
#include"../structures/domain/AbstractDomain.h"
#include"../structures/domain/SimpleDomain2D.h"
#include"../structures/domain/SimpleDomain.h"
#include"../structures/domain/IntersectionVascularizedDomain.h"
#include"../structures/domain/DomainNVR.h"
#include"../structures/domain/PartiallyVascularizedDomain.h"
#include"../structures/domain/DummyDomain.h"
#include"../structures/domain/StagedDomain.h"
#include"../constrains/AbstractConstraintFunction.h"
#include"../constrains/ConstantConstraintFunction.h"
#include"../constrains/ConstantPiecewiseConstraintFunction.h"
#include"../structures/tree/AbstractObjectCCOTree.h"
#include"../structures/tree/SingleVesselCCOOTree.h"
#include"../constrains/AbstractConstraintFunction.h"
#include"../core/MonitoringStagedFRROTreeGenerator.h"

class MonitoringStagedFRROTreeGeneratorLogger 
{
    // Attibutes
    MonitoringStagedFRROTreeGenerator *treeGenerator;
    FILE *file;

    public:
        /**
         * Constructor.
         * @param fp Pointer to file to write to, in "w" mode.
         * @param treeGen Tree generator.
         */
        MonitoringStagedFRROTreeGeneratorLogger(FILE *fp, MonitoringStagedFRROTreeGenerator* treeGen);
        ~MonitoringStagedFRROTreeGeneratorLogger();
        /**
         * Writes relevant parameters of @p treeGenerator_ to @p file.
         */
        void write();
        
    private:
        void logDomainFiles(FILE *fp, AbstractDomain *domain);
        void logCostEstimator(FILE *fp, AbstractCostEstimator *costEstimator);
        void logGenData(FILE *fp, GeneratorData *data);
        void logConstraint(FILE *fp, AbstractConstraintFunction<double, int> *constraint);
        void logConstraint(FILE *fp, AbstractConstraintFunction<double, double> *constraint);
        void logDomain(FILE *fp, AbstractDomain *domain, long long int n_term, AbstractConstraintFunction<double, int> *gam,
    AbstractConstraintFunction<double, int> *epsLim, AbstractConstraintFunction<double, int> *nu);
};

#endif // MONITORINGSTAGEDFRROTREEGENERATORLOGGER_H
