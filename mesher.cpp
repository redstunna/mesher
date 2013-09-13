#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <map>
#include <memory.h>

#include <Gmsh.h>
#include <GModel.h>
#include <MElement.h>
#include <MTetrahedron.h>
#include <MVertex.h>
#include <MFace.h>

#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkTetra.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkstd/string>

using namespace std;

// partitioning info consists of vector containing lists of
// tetrs for each zone in partition
typedef vector< vector<int> > partitionInfo;

typedef struct {
    // vertex index
    int num;
    // vetex coordinates
    float coords[3];
    // list of faces this vertex is member of
    vector<int> faces;
    // list of tetrs this vertex is member of
    vector<int> tetrs;
} vertex;

typedef struct {
    // face index
    int num;
    // face vertices
    int verts[3];
    // list of tetrs this vertex is member of
    vector<int> tetrs;
} face;

typedef struct {
    // tetr num
    int num;
    // list of faces this vertex is member of
    int faces[4];
} tetr;

typedef struct {
    // list of vertices
    vector<vertex> verts;
    // list of faces
    vector<face> faces;
    // list of tetrs
    vector<tetr> tetrs;
} grid;

// builds local grid representation basing on gmsh entity
grid buildGrid(GEntity* &ent) {
    partitionInfo pi;
    grid g;

    // iterate through volume elements and build our mesh representation
    set<int> verts_idx;
    map<int, int> new_idx;
    for (unsigned int i = 0; i < ent->getNumMeshElements(); i++) {
        MTetrahedron* t = (MTetrahedron*)ent->getMeshElement(i);
        // get list of vertices for tetr
        vector<MVertex*> v;
        t->getVertices(v);
        // iterate thru vertices and add them to grid
        for (int j = 0; j < 4; j++) {
            vertex nv;
            nv.num = v[j]->getNum();
            if (verts_idx.count(nv.num) == 0) {
                new_idx[nv.num] = g.verts.size();
                nv.coords[0] = v[j]->x();
                nv.coords[1] = v[j]->y();
                nv.coords[2] = v[j]->z();
                g.verts.push_back(nv);
                verts_idx.insert(nv.num);
            }
        }
        // iterate thru faces and add them to grid if they do not exist yet
        int fidx[4];
        for (int j = 0; j < 4; j++) {
            MFace f = t->getFace(j);
            vector<MVertex*> vt;
            f.getOrderedVertices(vt);
            bool found = false;
            for (unsigned int k = 0; k < g.faces.size(); k++) {
                if (g.faces[k].verts[0] == new_idx[vt[0]->getNum()] &&
                    g.faces[k].verts[1] == new_idx[vt[1]->getNum()] &&
                    g.faces[k].verts[2] == new_idx[vt[2]->getNum()]) {
                    // face already exist
                    fidx[j] = k;
                    found = true;
                    break;
                }

            }

            if (!found) {
                // face not found -- create new one
                face nf;
                nf.num = g.faces.size();
                fidx[j] = nf.num;
                nf.verts[0] = new_idx[vt[0]->getNum()];
                nf.verts[1] = new_idx[vt[1]->getNum()];
                nf.verts[2] = new_idx[vt[2]->getNum()];
                g.faces.push_back(nf);
            }
        }

        // add tetr to grid
        tetr nt;
        nt.num = i;
        memcpy(nt.faces, fidx, sizeof(int)*4);
        g.tetrs.push_back(nt);
    }

    cout << "Built local grid" << endl;
    cout << "Vertices: " << g.verts.size() << endl;
    cout << "Faces: " << g.faces.size() << endl;
    cout << "Tetrahedrons: " << g.tetrs.size() << endl;

    return g;
}

// performs grid partitioning
// this is very simple implementation that splits grid into `num` parts:
// parts 1..num-1 consist of one tetrahedron (k-th)
// part num consists of all rest tetrahedrons
//
// FIXME this function should be replaced with a proper implementation
partitionInfo partitionVolume(grid g, int num) {
    partitionInfo pi(num);

    for (int i = 0; i < num-1; i++)
        pi[i].push_back(i);
    for (int i = num-1; i < g.tetrs.size(); i++)
        pi[num-1].push_back(i);

    return pi;
}

// saves partitioned grid to bunch of vtk files
void dumpGrid(grid g, partitionInfo pi, string prefix) {
    // iterate thru partitions and save grids
    int pi_n = 0;
    for (auto &p: pi) {
        // create writer and grid
        vtkXMLUnstructuredGridWriter *gw = vtkXMLUnstructuredGridWriter::New();
        vtkUnstructuredGrid *ug = vtkUnstructuredGrid::New();
        // add points
        vtkPoints *pts = vtkPoints::New();
        for (auto &vert: g.verts)
            pts->InsertNextPoint(vert.coords[0], vert.coords[1], vert.coords[2]);
        ug->SetPoints(pts);
        // add cells
        vtkTetra *vt = vtkTetra::New();
        for (auto &idx: p) {
            set<int> idxs;
            for (int z = 0; z < 4; z++) {
                for (int zz = 0; zz < 3; zz++)
                    idxs.insert(g.faces[g.tetrs[idx].faces[z]].verts[zz]);
            }
            int n = 0;
            for (auto &k: idxs)
                vt->GetPointIds()->SetId(n++, k);
            ug->InsertNextCell(vt->GetCellType(), vt->GetPointIds());
        }

        // save dump
        gw->SetInput(ug);
        gw->SetFileName((prefix + "_" + to_string(pi_n++) + ".vtu").c_str());
        gw->Update();

        // clean up
        pts->Delete();
        vt->Delete();
        ug->Delete();
        gw->Delete();
    }
}

int main(int argc, char **argv) {
    // check if we have engough parameters
    if (argc != 3) {
        cerr << "Usage: mesher <geo file> <parts>" << endl;
        return -1;
    }

    // check if input file is readable
    ifstream f(argv[1]);
    if (!f.good()) {
        cerr << "Geo file does not exist." << endl;
        return -1;
    }

    // check if number of parts specified correctly
    int nparts;
    istringstream(argv[2]) >> nparts;
    if (nparts < 1) {
        cerr << "Invalid number of parts specified." << endl;
        return -1;
    }

    // init gmsh library
    GmshInitialize(1, argv);

    // set terminal output
    GmshSetOption("General", "Terminal", 1.0);
    // turn on mesh optimization
    GmshSetOption("Mesh", "Optimize", 1.0);

    // load model from geo file
    GModel m;
    m.setFactory("Gmsh");
    m.readGEO(argv[1]);

    // set meshing parameters basing on loaded object bounds
    auto bounds = m.bounds();
    GmshSetOption("Mesh", "CharacteristicLengthMin", bounds.diag()/8);
    GmshSetOption("Mesh", "CharacteristicLengthMax", bounds.diag()/4);

    // generate mesh
    m.mesh(3);

    // flush output stream
    cout.flush();
    cout << endl << endl << endl;

    // get mesh entities
    vector<GEntity*> ents;
    m.getEntities(ents);
    // iterate thru entities and partition volumes
    // FIXME add support for multiple volumes
    for (auto &ent: ents)
        if (ent->geomType() == GEntity::Volume) {
            // get partition information
            auto g = buildGrid(ent);
            auto pi = partitionVolume(g, nparts);
            // dump
            dumpGrid(g, pi, argv[1]);
            // stop generation (see FIXME above)
            break;
        }

    // finalize gmsh library
    GmshFinalize();

    return 0;
}
