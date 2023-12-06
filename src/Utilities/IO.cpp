#include "IO.h"

void IO::WriteVTU(Mesh &mesh, std::string filename) {
    Point *points = mesh.GetPoints();
    size_t pointCount = mesh.GetPointCount();
    size_t cellCount = mesh.GetCellCount();
    Cell *cellList = mesh.GetCellList();
    
    std::ofstream f;
    f.open(filename + ".vtu");
    f << "<VTKFile type=\"UnstructuredGrid\" version=\"2.2\" byte_order=\"LittleEndian\">\n";
    f << "  <UnstructuredGrid>\n"; 
    f << "    <Piece NumberOfPoints=" << "\"" << pointCount << "\"" << " NumberOfCells=" 
        << "\"" << cellCount << "\">\n"; 
    f << "      <Points>\n";
    f << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (size_t i = 0; i < pointCount; i++) {
        Point *p = &points[i];
        f << "          " << p->GetX() << " " << p->GetY() << " " << p->GetZ() << std::endl;
    }
    f << "        </DataArray>\n";
    f << "      </Points>\n";
    f << "      <Cells>\n";
    f << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (size_t i = 0; i < cellCount; i++) {
        Cell* cell = &cellList[i];
        f << "          " << cell->GetConnectivity() << "\n";
    }
    f << "        </DataArray>\n";
    f << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    int offset = 8;
    f << "          ";
    for (size_t i = 0; i < cellCount; i++) {
        f << offset << " ";
        offset += 8;
    }
    f << "\n        </DataArray>\n";
    f << "        <DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">\n";
    f << "          ";
    for (size_t i = 0; i < cellCount; i++) {
        f << 12 << " ";
    }
    f << "\n        </DataArray>\n";
    f << "      </Cells>\n";
    f << "      <CellData>\n";
    f << "        <DataArray type=\"Float64\" Name=\"U\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (size_t i = 0; i < cellCount; i++) {
        Cell* cell = &cellList[i];
        f << "          ";
        Vector3D vel = cell->GetVelocity();
        f << vel.GetU() << " " << vel.GetV() << " " << 0 << "\n";
    }
    f << "        </DataArray>\n";
    f << "        <DataArray type=\"Float64\" Name=\"p\" format=\"ascii\">\n";
    for (size_t i = 0; i < cellCount; i++) {
        Cell* cell = &cellList[i];
        f << "          ";
        double pressure = cell->GetPressure();
        f << pressure << "\n";
    }
    f << "        </DataArray>\n";
    f << "      </CellData>\n";
    f << "    </Piece>\n";
    f << "  </UnstructuredGrid>\n";
    f << "</VTKFile>\n";
    f.close();
}
