/// \sa ttk::ftc::FTCTree
#ifndef _VTK_CONTOURFORESTS_H
#define _VTK_CONTOURFORESTS_H

// VTK includes
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDataSetAlgorithm.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
// Data array
#include <vtkCharArray.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>

// Unused ? (compile without these on my computer)
#include <vtkFiltersCoreModule.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkLine.h>
#include <vtkObjectFactory.h>
#include <vtkType.h>
#include <vtkUnstructuredGrid.h>

// ttk code includes
#include <FTCTree.h>
#include <ttkWrapper.h>
#include<ttkFTCStructures.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkFTCTree : public vtkDataSetAlgorithm, public Wrapper
#else
class ttkFTCTree : public vtkDataSetAlgorithm, public Wrapper
#endif
{
  public:
   static ttkFTCTree* New();

   vtkTypeMacro(ttkFTCTree, vtkDataSetAlgorithm);

   // default ttk setters
   vtkSetMacro(debugLevel_, int);

   void SetThreadNumber(int threadNumber)
   {
      ThreadNumber = threadNumber;
      SetThreads();
   }

   void SetUseAllCores(bool onOff)
   {
      UseAllCores = onOff;
      SetThreads();
   }
   // end of default ttk setters

   vtkSetMacro(ScalarField, string);
   vtkGetMacro(ScalarField, string);

   vtkSetMacro(UseInputOffsetScalarField, int);
   vtkGetMacro(UseInputOffsetScalarField, int);

   vtkSetMacro(InputOffsetScalarFieldName, string);
   vtkGetMacro(InputOffsetScalarFieldName, string);

   vtkSetMacro(ScalarFieldId, int);
   vtkGetMacro(ScalarFieldId, int);

   vtkSetMacro(OffsetFieldId, int);
   vtkGetMacro(OffsetFieldId, int);

   // Parameters uses a structure, we can't use vtkMacro on them
   void SetTreeType(const int type)
   {
       params_.treeType = (ftc::TreeType)type;
       Modified();
   }

   ftc::TreeType GetTreeType(void) const
   {
       return params_.treeType;
   }

   void SetWithSegmentation(const bool segm)
   {
      params_.segm = segm;
      Modified();
   }

   bool GetWithSegmentation(void) const
   {
      return params_.segm;
   }

   void SetWithNormalize(const bool norm)
   {
      params_.normalize = norm;
      Modified();
   }

   bool GetWithNormalize(void) const
   {
      return params_.normalize;
   }

   void SetWithAdvStats(const bool adv)
   {
      params_.advStats = adv;
      Modified();
   }

   bool GetWithAdvStats(void) const
   {
      return params_.advStats;
   }

   void SetSuperArcSamplingLevel(int lvl)
   {
      params_.samplingLvl = lvl;
      Modified();
   }

   int GetSuperArcSamplingLevel(void) const
   {
      return params_.samplingLvl;
   }

   int setupTriangulation();
   int getScalars();
   int getOffsets();

   int getSkeletonNodes(vtkUnstructuredGrid* outputSkeletonNodes);

   int addDirectSkeletonArc(const ftc::idSuperArc arcId, const int cc, vtkPoints* points,
                            vtkUnstructuredGrid* skeletonArcs, ArcData& arcData);

   int addSampledSkeletonArc(const ftc::idSuperArc arcId, const int cc, vtkPoints* points,
                             vtkUnstructuredGrid* skeletonArcs, ArcData& arcData);

   int addCompleteSkeletonArc(const ftc::idSuperArc arcId, const int cc, vtkPoints* points,
                              vtkUnstructuredGrid* skeletonArcs, ArcData& arcData);

   int getSkeletonArcs(vtkUnstructuredGrid* outputSkeletonArcs);

   int getSegmentation(vtkDataSet* outputSegmentation);

#ifdef TTK_ENABLE_FTC_TREE_STATS_TIME
   void printCSVStats();
   void printCSVExtr(const ftc::FTCTree_MT* const tree) const;
   void printCSVTree(const ftc::FTCTree_MT* const tree, const bool printExtract = true) const;
#endif

  protected:
   ttkFTCTree();
   ~ttkFTCTree();

   TTK_SETUP();

   void identify(vtkDataSet* ds) const;

   virtual int FillInputPortInformation(int port, vtkInformation* info);
   virtual int FillOutputPortInformation(int port, vtkInformation* info);

  private:

   string ScalarField;
   bool   UseInputOffsetScalarField;
   string InputOffsetScalarFieldName;
   int    ScalarFieldId;
   int    OffsetFieldId;

   ftc::Params params_;

   int                            nbCC_;
   vector<vtkDataSet*>            connected_components_;
   vector<Triangulation*>         triangulation_;
   vector<LocalFTC>               ftcTree_;
   vector<vtkDataArray*>          inputScalars_;
   vector<vector<ftc::idVertex>>  offsets_;

   bool                   hasUpdatedMesh_;
};

#endif  // _VTK_CONTOURFORESTS_H
