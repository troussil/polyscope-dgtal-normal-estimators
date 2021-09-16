#include <iostream>
#include <chrono>
#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/helpers/Shortcuts.h>
#include <DGtal/helpers/ShortcutsGeometry.h>
#include <DGtal/topology/ImplicitDigitalSurface.h>
#include "polyscope/polyscope.h"
#include "polyscope/point_cloud.h"
#include "polyscope/surface_mesh.h"


namespace DGtal {
    namespace internal {
        using KSpace = Z3i::KSpace;
        using DigitizedImplicitShape = Shortcuts<KSpace>::DigitizedImplicitShape3D;
        using ImplicitSurfaceContainer = ImplicitDigitalSurface<KSpace, DigitizedImplicitShape>;
        using ImplicitDigitalSurface = DGtal::DigitalSurface<ImplicitSurfaceContainer>;
    } // namespace internal

    template < typename TPointPredicate >
    CountedPtr<internal::ImplicitDigitalSurface>
    makeImplicitDigitalSurface (CountedPtr<TPointPredicate> bimage,
                                internal::KSpace const& K,
                                Parameters const& params) {
        using SH3 = Shortcuts<internal::KSpace>;

        bool surfel_adjacency = params[ "surfelAdjacency" ].as<int>();
        int nb_tries_to_find_a_bel = params[ "nbTriesToFindABel" ].as<int>();

        CountedPtr<internal::ImplicitDigitalSurface> ptrSurface;
        SH3::Surfel bel;

        try {
            bel = Surfaces<internal::KSpace>::findABel(K, *bimage,
                                                       nb_tries_to_find_a_bel);
        } catch (DGtal::InputException& e) {
            return ptrSurface;
        }

        internal::ImplicitSurfaceContainer* container = new internal::ImplicitSurfaceContainer(K, *bimage,
                                                                                               surfel_adjacency, bel);

        return CountedPtr<internal::ImplicitDigitalSurface>(new internal::ImplicitDigitalSurface(container)); // acquired
    }
} // namespace DGtal


using namespace DGtal;
using namespace Z3i;

// Using standard 3D digital space.
typedef Shortcuts<Z3i::KSpace>         SH3;
typedef ShortcutsGeometry<Z3i::KSpace> SHG3;

typedef Shortcuts<KSpace> SH3;
typedef SH3::DigitalSurface Surface;
typedef SH3::Point Point;
typedef SH3::RealPoint RealPoint;
typedef SH3::KSpace KSpace;
typedef SH3::BinaryImage BinaryImage;
typedef Point::Coordinate Integer;

const Z3i::RealPoint shift = Z3i::RealPoint(-0.5, -0.5, -0.5); // for display
#include "PolyscopeUtils.hpp"

//global vars
float h=1.0;
float rradius=5.0;
const char* surfaces_equations[] = {
				    "3*x^2+2*y^2+z^2-90",
				    "x^2+y^2+z^2-100",
				    "(x^2+y^2+z^2+6*6-2*2)^2-4*6*6*(x^2+y^2)",
				    "2*x+6*y+15*z",
				    "x^2+y^2-z^2-90",
				    "(2*x+6*y+15*z)*(-3*x-4*y+3*z)",
				    "x^2+y^2+2*z^2-x*y*z+z^3-100",
				    "goursat",
				    "10000-(x^2+y^2+z^2+1000*(x^2+y^2)*(x^2+z^2)*(y^2+z^2))",
				    "100-(x^2*y^2*z^2+4*x^2+4*y^2+3*z^2)",
				    "x^2-(y^2+z^2)^2",
				    "-1*(x^2+2.25*y^2+z^2-1)^3+x^2*z^3+0.1125*y^2*z^3",
				    "-0.9*(y^2+z^2-1)^2-(x^2+y^2-1)^3"
};
const char* surfaces_names[] = {
				"ellipse",
				"sphere",
				"torus",
				"plane",
				"hyperboloid",
				"planes-inter",
				"sympa",
				"goursat",
				"distel",
				"leopold",
				"diabolo",
				"heart",
				"crixxi"
};
static int surfaces_current = 0;
static int img_current = 0; 
const std::string VOL_PATH = "../../data"; // TODO: relative to bin directory
const char* img_path[] = {
        "/rounded.vol",
        "/fandisk.vol",
        "/cubesphere.vol",
        "/torus.vol",
	"/rotated_cube.vol",
	"/cube.vol",
	"/INH5_64_0.vol",
	"/INH3_128_0.vol",
	"/E2bis_128_450-961_0.vol",
	"/INH5_128_bwliss_0.vol",
	"/top_sphere.vol",
	"/bunny-64.vol",	
	"/rotatedCube-2-6-15.vol"
};
const char* img_names[] = {
        "rounded",
        "fandisk",
        "cubesphere",
        "torus",
	"rotated cube",
	"cube", 
	"INH5_64_0",
	"INH3_128_0",
	"E2bis_128_450-961_0",
	"INH5_128_bwliss_0.vol",
	"top sphere",
	"bunny",
	"rotated cube 2 6 15"
};

// RealPoint getNormalizedVec(const Point& v) {
//   double n = v.norm(); 
//   return RealPoint( {v[0]/n, v[1]/n, v[2]/n} );
// }

void loadVol()
{
  std::cerr << VOL_PATH + img_path[img_current] << std::endl; 
  auto params = SH3::defaultParameters() | SHG3::defaultParameters() |  SHG3::parametersGeometryEstimation();
  auto binary_image= SH3::makeBinaryImage(VOL_PATH + img_path[img_current], params );
  auto K               = SH3::getKSpace( binary_image, params );
  auto surface         = SH3::makeLightDigitalSurface( binary_image, K, params );
  auto surfels         = SH3::getSurfelRange( surface, params );
  
  //Computing some differential quantities
  params("r-radius", rradius);
  params("t-ring", rradius);
  auto normals   = SHG3::getCTrivialNormalVectors(surface, surfels, params);
  auto normalsTrivial   = SHG3::getTrivialNormalVectors(K,surfels);
  auto normalsII = SHG3::getIINormalVectors(binary_image, surfels, params);
  auto Mcurv     = SHG3::getIIMeanCurvatures(binary_image, surfels, params);
  auto Gcurv     = SHG3::getIIGaussianCurvatures(binary_image, surfels, params);

  //Surfel area measure
  std::vector<double> areaMeasure(surfels.size());
  for(auto i=0; i < areaMeasure.size(); ++i)
    areaMeasure[i] = normalsTrivial[i].dot(normalsII[i]);

  // //
  // SH3::Cell2Index c2i;
  // auto primalSurface   = SH3::makePrimalPolygonalSurface(c2i, surface);
  
  // std::vector<std::vector<unsigned long>> faces;
  // for(auto &face: primalSurface->allFaces())
  //   faces.push_back(primalSurface->verticesAroundFace( face ));
  // auto digsurf = polyscope::registerSurfaceMesh("Primal", primalSurface->positions(), faces);
  
  // //Attaching quantities
  // digsurf->addFaceVectorQuantity("Trivial normal vectors", normalsTrivial);
  // digsurf->addFaceVectorQuantity("CTrivial normal vectors", normals);
  // digsurf->addFaceVectorQuantity("II normal vectors", normalsII);
  // digsurf->addFaceScalarQuantity("II mean curvature", Mcurv);
  // digsurf->addFaceScalarQuantity("II Gaussian curvature", Gcurv);
  // digsurf->addFaceScalarQuantity("Surfel area measure", areaMeasure);


  // Primal surface
  auto primalSurface = polyscopeSurfels(K, SH3::getSurfelRange(surface));
  auto digsurf = polyscope::registerSurfaceMesh("Primal surface",
						primalSurface.first,
						primalSurface.second);

  auto embedder        = SH3::getSCellEmbedder( K );
  std::vector<SH3::RealPoint> surfelCenters; 
  for (auto s: surfels) {
    surfelCenters.push_back(embedder(s)); 
  }
  auto pointCloud = polyscope::registerPointCloud("Surfel centers",
						  surfelCenters);
  
  //Attaching quantities
  pointCloud->addVectorQuantity("CTrivial normal vectors", normals);
  pointCloud->addVectorQuantity("II normal vectors", normalsII);


  
  // auto K           = SH3::getKSpace( binary_imageSnow );
  // auto surface    = SH3::makeLightDigitalSurfaces( binary_imageSnow, KSnow, params );
 
  // auto cpt=0;
  // trace.info()<<"Nb Surfaces= "<<surfacesSnow.size()<<std::endl;
  // for(auto &surf: surfacesSnow)
  //   {
  //     SH3::Cell2Index c2i;
  //     auto surfels         = SH3::getSurfelRange( surf, params );
  //     auto primalSurface   = SH3::makePrimalPolygonalSurface(c2i, surf);
  //     auto dualSurface = SH3::makeDualPolygonalSurface(surf);
    
    
  //     if (surfels.size()< 10) continue;
  //     //    trace.info()<<"Nb surfels =" <<surfels.size()<<std::endl;
  //     //    trace.info()<<"Nb Faces =" <<dualSurface->allFaces().size()<<std::endl;
  //     //    trace.info()<<"K " <<KSnow<<std::endl;
  //     //
  //     //Need to convert the faces
  //     std::vector<std::vector<unsigned long>> faces;
  //     for(auto &face: dualSurface->allFaces())
  // 	faces.push_back(dualSurface->verticesAroundFace( face ));
    
  //     auto name ="Dual surface "+std::to_string(cpt);
  //     cpt++;
    
  //     polyscope::registerSurfaceMesh(name, dualSurface->positions(), faces);
  //   }
}


// void otherQuantities()
// {

//   auto params = SH3::defaultParameters() | SHG3::defaultParameters() |  SHG3::parametersGeometryEstimation();
//   params( "polynomial", "goursat" )( "gridstep", h );
//   auto implicit_shape  = SH3::makeImplicitShape3D  ( params );
//   auto digitized_shape = SH3::makeDigitizedImplicitShape3D( implicit_shape, params );
//   auto binary_image    = SH3::makeBinaryImage( digitized_shape, params );
//   auto K               = SH3::getKSpace( params );
//   auto embedder        = SH3::getCellEmbedder( K );

//   CountedPtr<Surface> surface = SH3::makeDigitalSurface(binary_image, K, params);

//   SH3::Cell2Index c2i;
//   auto primalSurface   = SH3::makePrimalPolygonalSurface(c2i, surface);
  
//   std::vector<std::vector<unsigned long>> primal_faces;
//   for(auto &face: primalSurface->allFaces())
//     primal_faces.push_back(primalSurface->verticesAroundFace( face ));
//   auto digsurf = polyscope::registerSurfaceMesh("Primal", primalSurface->positions(), primal_faces);
//   //digsurf->edgeWidth=1.0;
//   //digsurf->edgeColor={1.,1.,1.};

//   // Iterate over the surfels
//   auto surfels = SH3::getSurfelRange(surface, params);
//   std::clog << surfels.size() << " surfels" << std::endl;
  
//   // //Attaching quantities
//   // digsurf->addFaceVectorQuantity("H normal vectors", normals);

//   // //Attaching polygon soup
//   // std::vector<RealPoint> positions;
//   // std::vector<std::vector<unsigned long> > faces;
//   // auto approx = polyscope::registerSurfaceMesh("Approximation", positions, faces);
  
// }

void updateQuantities()
{
  auto params = SH3::defaultParameters() | SHG3::defaultParameters() |  SHG3::parametersGeometryEstimation();
  //params( "polynomial", "goursat" )( "gridstep", h );
  params( "polynomial", surfaces_equations[surfaces_current] )( "gridstep", h );

  // // Build surface
  // std::chrono::system_clock::time_point before_dis = std::chrono::system_clock::now();
  // auto implicit_shape  = SH3::makeImplicitShape3D  ( params );
  // auto digitized_shape = SH3::makeDigitizedImplicitShape3D( implicit_shape, params );
  // std::chrono::system_clock::time_point after_dis = std::chrono::system_clock::now();
  // std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(after_dis - before_dis).count() << " ms" << std::endl;
  // KSpace K = SH3::getKSpace(params);
  // std::chrono::system_clock::time_point before_surf = std::chrono::system_clock::now();
  // auto surface = DGtal::makeImplicitDigitalSurface(digitized_shape, K, params);
  // std::chrono::system_clock::time_point after_surf = std::chrono::system_clock::now();
  // std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(after_surf - before_surf).count() << " ms" << std::endl;

  // auto surfels = SH3::getSurfelRange(surface);
  
  auto implicit_shape  = SH3::makeImplicitShape3D  ( params );
  auto digitized_shape = SH3::makeDigitizedImplicitShape3D( implicit_shape, params );
  auto binary_image    = SH3::makeBinaryImage( digitized_shape, params );
  auto K               = SH3::getKSpace( params );
  auto surface         = SH3::makeLightDigitalSurface( binary_image, K, params );
  auto surfels         = SH3::getSurfelRange( surface, params );
  
  //Computing some differential quantities
  params("r-radius", rradius);
  auto normalsTruth = SHG3::getNormalVectors(implicit_shape, K, surfels, params);
  auto normals   = SHG3::getCTrivialNormalVectors(surface, surfels, params);
  auto normalsTrivial   = SHG3::getTrivialNormalVectors(K,surfels);
  auto normalsII = SHG3::getIINormalVectors(binary_image, surfels, params);
  auto Mcurv     = SHG3::getIIMeanCurvatures(binary_image, surfels, params);
  auto Gcurv     = SHG3::getIIGaussianCurvatures(binary_image, surfels, params);

  //Surfel area measure
  std::vector<double> areaMeasure(surfels.size());
  for(auto i=0; i < areaMeasure.size(); ++i)
    areaMeasure[i] = normalsTrivial[i].dot(normalsII[i]);

  //Normal error
  std::vector<double> normalsError(surfels.size());
  for(auto i=0; i < normalsError.size(); ++i)
    normalsError[i] = normalsTruth[i].dot(normalsII[i]);

  //
  SH3::Cell2Index c2i;
  auto primalSurface   = SH3::makePrimalPolygonalSurface(c2i, surface);
  
  std::vector<std::vector<unsigned long>> faces;
  for(auto &face: primalSurface->allFaces())
    faces.push_back(primalSurface->verticesAroundFace( face ));
  auto digsurf = polyscope::registerSurfaceMesh("Primal", primalSurface->positions(), faces);

  //Attaching quantities
  digsurf->addFaceVectorQuantity("Truth normals", normalsTruth);
  digsurf->addFaceVectorQuantity("Trivial normal vectors", normalsTrivial);
  digsurf->addFaceVectorQuantity("CTrivial normal vectors", normals);
  digsurf->addFaceVectorQuantity("II normal vectors", normalsII);
  digsurf->addFaceScalarQuantity("II mean curvature", Mcurv);
  digsurf->addFaceScalarQuantity("II Gaussian curvature", Gcurv);
  digsurf->addFaceScalarQuantity("Surfel area measure", areaMeasure);
  digsurf->addFaceScalarQuantity("II normal error", normalsError);

}

// void updateGeom()
// {
//   auto params = SH3::defaultParameters() | SHG3::defaultParameters() |  SHG3::parametersGeometryEstimation();
  
//   params("surfaceComponents", "All");
//   //params( "polynomial", "goursat" )( "gridstep", h );
//   params( "polynomial", surfaces_equations[surfaces_current] )( "gridstep", h );
//   auto implicit_shape  = SH3::makeImplicitShape3D  ( params );
//   auto digitized_shape = SH3::makeDigitizedImplicitShape3D( implicit_shape, params );
//   auto binary_image    = SH3::makeBinaryImage( digitized_shape, params );
//   auto K               = SH3::getKSpace( params );
//   auto embedder        = SH3::getCellEmbedder( K );
//   auto surface         = SH3::makeLightDigitalSurface( binary_image, K, params );
  
//   //Polynomial suface
//   //Need to convert the faces
//   SH3::Cell2Index c2i;
//   auto surfels         = SH3::getSurfelRange( surface, params );
//   auto primalSurface   = SH3::makePrimalPolygonalSurface(c2i, surface);
  
//   std::vector<std::vector<unsigned long>> faces;
//   for(auto &face: primalSurface->allFaces())
//     faces.push_back(primalSurface->verticesAroundFace( face ));
//   auto digsurf = polyscope::registerSurfaceMesh("Primal", primalSurface->positions(), faces);
// }

// Your callback functions
void myCallback()
{
  ImGui::Combo("surface", &surfaces_current,
	       surfaces_names, IM_ARRAYSIZE(surfaces_names));

  
  ImGui::SliderFloat("h", &h, 0.001, 1);

  ImGui::Combo("volume", &img_current,
	       img_names, IM_ARRAYSIZE(img_names));

  ImGui::SliderFloat("r-radius", &rradius, 0, 10);
  
  if (ImGui::Button("Update"))
    updateQuantities();
  // if (ImGui::Button("???"))
  //   otherQuantities();
  if (ImGui::Button("Load 3d volume"))
    loadVol();
}

int main()
{
  polyscope::init();
  
  // Specify the callback
  polyscope::state::userCallback = myCallback;
  polyscope::show();
  return 0;
}
