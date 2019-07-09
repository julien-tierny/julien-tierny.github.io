/*
 * file:                  Editor.cpp
 * description:           Data-structures and processing.
 * author:                Your Name Here <Your Email Address Here>.
 * date:                  The Date Here.
 */

#include "Editor.h"

#include <vtksys/SystemTools.hxx>

Editor::Editor()
{
   grid_                  = NULL;
   lastObject_            = true;
   debug_                 = 1;
   ftcTree_               = ttkFTCTree::New();
   core_                  = 1 ;
   fieldId_               = 0;
   treeType_              = 2;
   ttk::globalDebugLevel_ = 3;
}

Editor::~Editor()
{
   // delete the mesh reader and the associated data at once
   ftcTree_->Delete();
}

int Editor::execute()
{
   ftcTree_->setDebugLevel(debug_);
   ftcTree_->SetdebugLevel_(debug_);
   ftcTree_->SetUseAllCores(false);
   ftcTree_->SetThreadNumber(core_);

   ftcTree_->SetInputData(grid_);
   ftcTree_->SetScalarFieldId(fieldId_);
   ftcTree_->SetTreeType(treeType_);

   ftcTree_->Update();

   grid_->ShallowCopy(ftcTree_->GetOutput());

   return 0;
}

int Editor::init(int &argc, char **argv)
{
   CommandLineParser parser;

   // specify argument "-g" : The path of the grid: REQUIRED
   parser.setArgument("g", &inputFilePath_, "Path to the input 3D grid");
   parser.setArgument("f", &fieldId_, "Field identifier", true);
   parser.setArgument("T", &treeType_, "type of tree : 2 is CT", true);

   // now parse the command line
   parser.parse(argc, argv);

   // set default values
   debug_ = ttk::globalDebugLevel_;
   core_ = parser.getThreadNumber();

   // now load the data to the editor
   loadData();

   return 0;
}

int Editor::loadData()
{
   // create a reader object
   std::string extension = vtksys::SystemTools::GetFilenameLastExtension(inputFilePath_);

   if (extension == ".vtu") {
      grid_ = ReadAnXMLFile<vtkXMLUnstructuredGridReader>(inputFilePath_.c_str());
   } else if (extension == ".vti") {
      grid_ = ReadAnXMLFile<vtkXMLImageDataReader>(inputFilePath_.c_str());
   } else {
      cerr << "Bad format, need vtu" << endl;
      return -1;
   }

   {
      stringstream msg;
      msg << "[Editor]   done! (read " << grid_->GetNumberOfPoints() << " vertices, "
          << grid_->GetNumberOfCells() << " cells)" << endl;
      dMsg(cout, msg.str(), 1);
   }

   return 0;
}

int Editor::saveData(const string &fileName) const
{
   vtkXMLUnstructuredGridWriter *writer = vtkXMLUnstructuredGridWriter::New();

   writer->SetFileName(fileName.data());
   writer->SetInputData((vtkUnstructuredGrid *)ftcTree_->GetOutput());
   writer->Write();

   writer->Delete();

   return 0;
}
