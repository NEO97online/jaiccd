/***
 * libccd
 * ---------------------------------
 * Copyright (c)2010 Daniel Fiser <danfis@danfis.cz>
 *
 *
 *  This file is part of libccd.
 *
 *  Distributed under the OSI-approved BSD License (the "License");
 *  see accompanying file BDS-LICENSE for details or see
 *  <http://www.opensource.org/licenses/bsd-license.php>.
 *
 *  This software is distributed WITHOUT ANY WARRANTY; without even the
 *  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *  See the License for more information.
 */

CCD_PT_VERTEX :: 1;
CCD_PT_EDGE   :: 2;
CCD_PT_FACE   :: 3;

/**
 * General polytope element.
 * Could be vertex, edge or triangle.
 */
CCD_PolytopeElement :: struct {
    type: int;           /*! type of element */
    dist: ccd_real_t;    /*! distance from origin */
    witness: Vector3;    /*! witness point of projection of origin */
    list: CCD_List;    /*! list of elements of same type */
}

/**
 * Polytope's vertex.
 */
CCD_PolytopeVertex :: struct {
    using base: CCD_PolytopeElement;

    id: int;
    v: CCD_Support;
    edges: CCD_List; //!< List of edges
}

/**
 * Polytope's edge.
 */
CCD_PolytopeEdge :: struct {
    using base: CCD_PolytopeElement;

    vertex: [2] *CCD_PolytopeVertex; //!< Reference to vertices
    faces: [2] *CCD_PolytopeFace; //!< Reference to faces

    vertex_list: [2] CCD_List; //!< List items in vertices' lists
}

/**
 * Polytope's triangle faces.
 */
CCD_PolytopeFace :: struct {
    using base: CCD_PolytopeElement;

    edge: [3] *CCD_PolytopeEdge; //!< Reference to surrounding edges
}

/**
 * Struct containing polytope.
 */
CCD_Polytope :: struct {
    vertices: CCD_List; //!< List of vertices
    edges: CCD_List; //!< List of edges
    faces: CCD_List; //!< List of faces

    nearest: *CCD_PolytopeElement;
    nearest_dist: ccd_real_t;
    nearest_type: int;
}

/**** INLINES ****/
/**
 * Adds vertex to polytope and returns pointer to newly created vertex.
 */
ccdPtAddVertexCoords :: inline (pt: *CCD_Polytope, x: ccd_real_t, y: ccd_real_t, z: ccd_real_t) -> *CCD_PolytopeVertex {
    s: CCD_Support;
    ccdVec3Set(*s.v, x, y, z);
    return ccdPtAddVertex(pt, *s);
}

/**
 * Deletes vertex from polytope.
 * Returns 0 on success, -1 otherwise.
 */
ccdPtDelVertex :: inline (pt: *CCD_Polytope, v: *CCD_PolytopeVertex) -> int {
    // test if any edge is connected to this vertex
    if (!ccdListEmpty(*v.edges)) {
        assert(false);
        return -1;
    }

    assert(v.type == CCD_PT_VERTEX);

    // delete vertex from main list
    ccdListDel(*v.list);

    if (cast(*void)(pt.nearest) == cast(*void)v){
        pt.nearest = null;
    }

    verts_removed += 1;

    // free(v);
    return 0;
}

ccdPtDelEdge :: inline (pt: *CCD_Polytope, e: *CCD_PolytopeEdge) -> int {
    // text if any face is connected to this edge (faces[] is always
    // aligned to lower indices)
    if (e.faces[0] != null) {
        assert(false);
        return -1;
    }

    assert(e.type == CCD_PT_EDGE);

    // disconnect edge from lists of edges in vertex struct
    ccdListDel(*e.vertex_list[0]);
    ccdListDel(*e.vertex_list[1]);

    // disconnect edge from main list
    ccdListDel(*e.list);

    if (cast(*void)(pt.nearest) == cast(*void)e){
        pt.nearest = null;
    }

    edges_removed += 1;

    // free(e);
    return 0;
}

ccdPtDelFace :: inline (pt: *CCD_Polytope, f: *CCD_PolytopeFace) -> int {
    e: *CCD_PolytopeEdge;
    i: size_t;

    assert(f.type == CCD_PT_FACE);

    // remove face from edges' recerence lists
    for i : 0..2 { //(i = 0; i < 3; i++)
        e = f.edge[i];
        if (e.faces[0] == f){
            e.faces[0] = e.faces[1];
        }
        e.faces[1] = null;
    }

    // remove face from list of all faces
    ccdListDel(*f.list);

    if (cast(*void)(pt.nearest) == cast(*void)f){
        pt.nearest = null;
    }

    faces_removed += 1;

    // free(f);
    return 0;
}

/**
 * Returns vertices surrounding given triangle face.
 */
ccdPtFaceVec3 :: inline (face: *CCD_PolytopeFace, a: **Vector3, b: **Vector3, c: **Vector3) {
    a.* = *face.edge[0].vertex[0].v.v;
    b.* = *face.edge[0].vertex[1].v.v;

    if (face.edge[1].vertex[0] != face.edge[0].vertex[0]
            && face.edge[1].vertex[0] != face.edge[0].vertex[1]){
        c.* = *face.edge[1].vertex[0].v.v;
    }else{
        c.* = *face.edge[1].vertex[1].v.v;
    }
}

ccdPtFaceVertices :: inline (face: *CCD_PolytopeFace, a: **CCD_PolytopeVertex, b: **CCD_PolytopeVertex, c: **CCD_PolytopeVertex) {
    a.* = face.edge[0].vertex[0];
    b.* = face.edge[0].vertex[1];

    if (face.edge[1].vertex[0] != face.edge[0].vertex[0]
            && face.edge[1].vertex[0] != face.edge[0].vertex[1]){
        c.* = face.edge[1].vertex[0];
    }else{
        c.* = face.edge[1].vertex[1];
    }
}

ccdPtFaceEdges :: inline (f: *CCD_PolytopeFace, a: **CCD_PolytopeEdge, b: **CCD_PolytopeEdge, c: **CCD_PolytopeEdge) {
    a.* = f.edge[0];
    b.* = f.edge[1];
    c.* = f.edge[2];
}

ccdPtEdgeVec3 :: inline (e: *CCD_PolytopeEdge, a: **Vector3, b: **Vector3) {
    a.* = *e.vertex[0].v.v;
    b.* = *e.vertex[1].v.v;
}

ccdPtEdgeVertices :: inline (e: *CCD_PolytopeEdge, a: **CCD_PolytopeVertex, b: **CCD_PolytopeVertex) {
    a.* = e.vertex[0];
    b.* = e.vertex[1];
}

ccdPtEdgeFaces :: inline (e: *CCD_PolytopeEdge, f1: **CCD_PolytopeFace, f2: **CCD_PolytopeFace) {
    f1.* = e.faces[0];
    f2.* = e.faces[1];
}

_ccdPtNearestUpdate :: inline (pt: *CCD_Polytope, el: *CCD_PolytopeElement) {
    nearest := pt.nearest_dist;
    el_dist := el.dist;
    if (ccdEq(nearest, el_dist)){
        if (el.type < pt.nearest_type){
            pt.nearest = el;
            pt.nearest_dist = el.dist;
            pt.nearest_type = el.type;
        }
    }else if (el.dist < pt.nearest_dist){
        pt.nearest = el;
        pt.nearest_dist = el.dist;
        pt.nearest_type = el.type;
    }
}

_ccdPtNearestRenew :: (pt: *CCD_Polytope) {
    v: *CCD_PolytopeVertex;
    e: *CCD_PolytopeEdge;
    f: *CCD_PolytopeFace;

    pt.nearest_dist = CCD_REAL_MAX;
    pt.nearest_type = 3;
    pt.nearest = null;

    ccdListForEachEntry(*pt.vertices, *v, CCD_PolytopeVertex, #code {
        _ccdPtNearestUpdate(pt, cast(*CCD_PolytopeElement)v);
    });

    ccdListForEachEntry(*pt.edges, *e, CCD_PolytopeEdge, #code {
        _ccdPtNearestUpdate(pt, cast(*CCD_PolytopeElement)e);
    });

    ccdListForEachEntry(*pt.faces, *f, CCD_PolytopeFace, #code {
        _ccdPtNearestUpdate(pt, cast(*CCD_PolytopeElement)f);
    });
}

ccdPtInit :: (pt: *CCD_Polytope) {
    ccdListInit(*pt.vertices);
    ccdListInit(*pt.edges);
    ccdListInit(*pt.faces);

    pt.nearest = null;
    pt.nearest_dist = CCD_REAL_MAX;
    pt.nearest_type = 3;
}

ccdPtDestroy :: (pt: *CCD_Polytope) {
    f, f2: *CCD_PolytopeFace;
    e, e2: *CCD_PolytopeEdge;
    v, v2: *CCD_PolytopeVertex;

    // first delete all faces
    f = ccdListEntry(pt.faces.next, CCD_PolytopeFace);
    f2 = ccdListEntry(f.list.next, CCD_PolytopeFace);
    while *f.list != *pt.faces {
        ccdPtDelFace(pt, f);
        f = f2;
        f2 = ccdListEntry(f2.list.next, CCD_PolytopeFace);
    }

    // delete all edges
    e = ccdListEntry(pt.edges.next, CCD_PolytopeEdge);
    e2 = ccdListEntry(e.list.next, CCD_PolytopeEdge);
    while *e.list != *pt.edges {
        ccdPtDelEdge(pt, e);
        e = e2;
        e2 = ccdListEntry(e2.list.next, CCD_PolytopeEdge);
    }

    // delete all vertices
    v = ccdListEntry(pt.vertices.next, CCD_PolytopeVertex);
    v2 = ccdListEntry(v.list.next, CCD_PolytopeVertex);
    while *v.list != *pt.vertices {
        ccdPtDelVertex(pt, v);
        v = v2;
        v2 = ccdListEntry(v2.list.next, CCD_PolytopeVertex);
    }
}


ccdPtAddVertex :: (pt: *CCD_Polytope, v: *CCD_Support) -> *CCD_PolytopeVertex {
    vert: *CCD_PolytopeVertex;

    vert = New(CCD_PolytopeVertex, allocator=temp);
    if (vert == null)
        return null;

    vert.type = CCD_PT_VERTEX;
    ccdSupportCopy(*vert.v, v);

    vert.dist = ccdVec3Len2(*vert.v.v);
    ccdVec3Copy(*vert.witness, *vert.v.v);

    ccdListInit(*vert.edges);

    // add vertex to list
    ccdListAppend(*pt.vertices, *vert.list);

    // update position in .nearest array
    _ccdPtNearestUpdate(pt, cast(*CCD_PolytopeElement)vert);

    verts_added += 1;

    return vert;
}

/**
 * Adds edge to polytope.
 */
ccdPtAddEdge :: (pt: *CCD_Polytope, v1: *CCD_PolytopeVertex, v2: *CCD_PolytopeVertex) -> *CCD_PolytopeEdge {
    a, b: *Vector3;
    edge: *CCD_PolytopeEdge;

    if (v1 == null || v2 == null)
        return null;

    edge = New(CCD_PolytopeEdge, allocator=temp);
    if (edge == null)
        return null;

    edge.type = CCD_PT_EDGE;
    edge.vertex[0] = v1;
    edge.vertex[1] = v2;
    edge.faces[0], edge.faces[1] = null;

    a = *edge.vertex[0].v.v;
    b = *edge.vertex[1].v.v;
    edge.dist = ccdVec3PointSegmentDist2(*ccd_vec3_origin, a, b, *edge.witness);

    ccdListAppend(*edge.vertex[0].edges, *edge.vertex_list[0]);
    ccdListAppend(*edge.vertex[1].edges, *edge.vertex_list[1]);

    ccdListAppend(*pt.edges, *edge.list);

    // update position in .nearest array
    _ccdPtNearestUpdate(pt, cast(*CCD_PolytopeElement)edge);

    edges_added += 1;

    return edge;
}

/**
 * Adds face to polytope.
 */
ccdPtAddFace :: (pt: *CCD_Polytope, e1: *CCD_PolytopeEdge, e2: *CCD_PolytopeEdge, e3: *CCD_PolytopeEdge) -> *CCD_PolytopeFace {
    a, b, c: *Vector3;
    face: *CCD_PolytopeFace;
    e: *CCD_PolytopeEdge;
    i: size_t;

    if (e1 == null || e2 == null || e3 == null)
        return null;

    face = New(CCD_PolytopeFace, allocator=temp);
    if (face == null)
        return null;

    face.type = CCD_PT_FACE;
    face.edge[0] = e1;
    face.edge[1] = e2;
    face.edge[2] = e3;

    // obtain triplet of vertices
    a = *face.edge[0].vertex[0].v.v;
    b = *face.edge[0].vertex[1].v.v;
    e = face.edge[1];
    if (e.vertex[0] != face.edge[0].vertex[0]
            && e.vertex[0] != face.edge[0].vertex[1]){
        c = *e.vertex[0].v.v;
    }else{
        c = *e.vertex[1].v.v;
    }
    face.dist = ccdVec3PointTriDist2(*ccd_vec3_origin, a, b, c, *face.witness);

    for i : 0..2 { // (i = 0; i < 3; i++){
        if (face.edge[i].faces[0] == null){
            face.edge[i].faces[0] = face;
        }else{
            face.edge[i].faces[1] = face;
        }
    }

    ccdListAppend(*pt.faces, *face.list);

    // update position in .nearest array
    _ccdPtNearestUpdate(pt, cast(*CCD_PolytopeElement)face);

    faces_added += 1;

    return face;
}

/**
 * Recompute distances from origin for all elements in pt.
 */
ccdPtRecomputeDistances :: (pt: *CCD_Polytope) {
    v: *CCD_PolytopeVertex;
    e: *CCD_PolytopeEdge;
    f: *CCD_PolytopeFace;
    a, b, c: *Vector3;
    dist: ccd_real_t;

    ccdListForEachEntry(*pt.vertices, v, CCD_PolytopeVertex, list, #code {
        dist = ccdVec3Len2(*v.v.v);
        v.dist = dist;
        ccdVec3Copy(*v.witness, *v.v.v);
    });

    ccdListForEachEntry(*pt.edges, e, CCD_PolytopeEdge, list, #code {
        a = *e.vertex[0].v.v;
        b = *e.vertex[1].v.v;
        dist = ccdVec3PointSegmentDist2(ccd_vec3_origin, a, b, *e.witness);
        e.dist = dist;
    });

    ccdListForEachEntry(*pt.faces, f, CCD_PolytopeFace, list, #code {
        // obtain triplet of vertices
        a = *f.edge[0].vertex[0].v.v;
        b = *f.edge[0].vertex[1].v.v;
        e = f.edge[1];
        if (e.vertex[0] != f.edge[0].vertex[0]
                && e.vertex[0] != f.edge[0].vertex[1]){
            c = *e.vertex[0].v.v;
        }else{
            c = *e.vertex[1].v.v;
        }

        dist = ccdVec3PointTriDist2(ccd_vec3_origin, a, b, c, *f.witness);
        f.dist = dist;
    });
}

/**
 * Returns nearest element to origin.
 */
ccdPtNearest :: (pt: *CCD_Polytope) -> *CCD_PolytopeElement {
    if (!pt.nearest){
        _ccdPtNearestRenew(pt);
    }
    return pt.nearest;
}


// JAI: these are just debug functions i think so we probably dont need em rn.

// ccdPtDumpSVT :: (pt: *CCD_Polytope, fn: *char) { // fn should be what?? it was a `const char *fn`
//     FILE *fout;

//     fout = fopen(fn, "a");
//     if (fout == null)
//         return;

//     ccdPtDumpSVT2(pt, fout);

//     fclose(fout);
// }

// void ccdPtDumpSVT2(CCD_Polytope *pt, FILE *fout)
// {
//     CCD_PolytopeVertex *v, *a, *b, *c;
//     CCD_PolytopeEdge *e;
//     CCD_PolytopeFace *f;
//     size_t i;

//     fprintf(fout, "-----\n");

//     fprintf(fout, "Points:\n");
//     i = 0;
//     ccdListForEachEntry(&pt->vertices, v, CCD_PolytopeVertex, list){
//         v->id = i++;
//         fprintf(fout, "%lf %lf %lf\n",
//                 ccdVec3X(&v->v.v), ccdVec3Y(&v->v.v), ccdVec3Z(&v->v.v));
//     }

//     fprintf(fout, "Edges:\n");
//     ccdListForEachEntry(&pt->edges, e, CCD_PolytopeEdge, list){
//         fprintf(fout, "%d %d\n", e->vertex[0]->id, e->vertex[1]->id);
//     }

//     fprintf(fout, "Faces:\n");
//     ccdListForEachEntry(&pt->faces, f, CCD_PolytopeFace, list){
//         a = f->edge[0]->vertex[0];
//         b = f->edge[0]->vertex[1];
//         c = f->edge[1]->vertex[0];
//         if (c == a || c == b){
//             c = f->edge[1]->vertex[1];
//         }
//         fprintf(fout, "%d %d %d\n", a->id, b->id, c->id);
//     }
// }