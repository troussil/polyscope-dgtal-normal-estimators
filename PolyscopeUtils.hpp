#pragma once

// For the separating polygons, convert it to polyscope style
std::pair<std::vector<RealPoint>, std::vector<std::vector<unsigned long>>>
polygonsToPolyscope (std::vector<std::vector<RealPoint>> const& polygons) {
    std::vector<RealPoint> vertices;
    std::vector<std::vector<unsigned long>> faces;

    unsigned long idx_vertex = 0;
    for (const auto& polygon: polygons) {
        std::vector<unsigned long> face;
        for (const auto& v: polygon) {
            vertices.push_back(v + shift);
            face.push_back(idx_vertex++);
        }
        faces.push_back(face);
    }

    return { vertices, faces };
}

// Some meshes are not manifold (intersection of two planes for instance), we draw each surfel
// as a quad.
std::pair<std::vector<RealPoint>, std::vector<std::vector<unsigned long>>>
polyscopeSurfels (KSpace const& K, std::vector<SH3::Surfel> const& surfels) {
    std::vector<RealPoint> vertices;
    std::vector<std::vector<unsigned long>> faces;

    unsigned long idx_vertex = 0;
    for (const auto& s: surfels) {
        auto surfel_vertices = SH3::getPrimalVertices(K, s, true);

        vertices.push_back(K.uCoords(surfel_vertices[0]) + shift); // idx_vertex
        vertices.push_back(K.uCoords(surfel_vertices[1]) + shift); // idx_vertex + 1
        vertices.push_back(K.uCoords(surfel_vertices[2]) + shift); // idx_vertex + 2
        vertices.push_back(K.uCoords(surfel_vertices[3]) + shift); // idx_vertex + 3

        faces.push_back( { idx_vertex, idx_vertex + 1, idx_vertex + 2, idx_vertex + 3 });

        idx_vertex += 4;
    }

    return { vertices, faces };
}

