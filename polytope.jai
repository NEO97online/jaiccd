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


__CCD_PT_EL :: struct {
    type: int;           /*! type of element */
    dist: ccd_real_t;    /*! distance from origin */
    witness: Vector3;    /*! witness point of projection of origin */
    list: ccd_list_t;    /*! list of elements of same type */
}

/**
 * General polytope element.
 * Could be vertex, edge or triangle.
 */
ccd_pt_el_t :: struct {
    using base: __CCD_PT_EL;
}

/**
 * Polytope's vertex.
 */
ccd_pt_vertex_t :: struct {
    using base: __CCD_PT_EL;

    id: int;
    v: ccd_support_t;
    edges: ccd_list_t; //!< List of edges
}

/**
 * Polytope's edge.
 */
ccd_pt_edge_t :: struct {
    using base: __CCD_PT_EL;

    vertex: [2] *ccd_pt_vertex_t; //!< Reference to vertices
    faces: [2] *ccd_pt_face_t; //!< Reference to faces

    vertex_list: [2] ccd_list_t; //!< List items in vertices' lists
}

/**
 * Polytope's triangle faces.
 */
ccd_pt_face_t :: struct {
    using base: __CCD_PT_EL;

    edge: [3] *ccd_pt_edge_t; //!< Reference to surrounding edges
}

/**
 * Struct containing polytope.
 */
ccd_pt_t :: struct {
    vertices: ccd_list_t; //!< List of vertices
    edges: ccd_list_t; //!< List of edges
    faces: ccd_list_t; //!< List of faces

    nearest: *ccd_pt_el_t;
    nearest_dist: ccd_real_t;
    nearest_type: int;
}

/**** INLINES ****/
/**
 * Adds vertex to polytope and returns pointer to newly created vertex.
 */
ccdPtAddVertexCoords :: inline (pt: *ccd_pt_t, x: ccd_real_t, y: ccd_real_t, z: ccd_real_t) -> *ccd_pt_vertex_t {
    s: ccd_support_t;
    ccdVec3Set(*s.v, x, y, z);
    return ccdPtAddVertex(pt, *s);
}

/**
 * Deletes vertex from polytope.
 * Returns 0 on success, -1 otherwise.
 */
ccdPtDelVertex :: inline (pt: *ccd_pt_t, v: *ccd_pt_vertex_t) -> int {
    // test if any edge is connected to this vertex
    if (!ccdListEmpty(*v.edges))
        return -1;

    // delete vertex from main list
    ccdListDel(*v.list);

    if (cast(*void)(pt.nearest) == cast(*void)v){
        pt.nearest = null;
    }

    free(v);
    return 0;
}

ccdPtDelEdge :: inline (pt: *ccd_pt_t, e: *ccd_pt_edge_t) -> int {
    // text if any face is connected to this edge (faces[] is always
    // aligned to lower indices)
    if (e.faces[0] != null)
        return -1;

    // disconnect edge from lists of edges in vertex struct
    ccdListDel(*e.vertex_list[0]);
    ccdListDel(*e.vertex_list[1]);

    // disconnect edge from main list
    ccdListDel(*e.list);

    if (cast(*void)(pt.nearest) == cast(*void)e){
        pt.nearest = null;
    }

    free(e);
    return 0;
}

ccdPtDelFace :: inline (pt: *ccd_pt_t, f: *ccd_pt_face_t) -> int {
    e: *ccd_pt_edge_t;
    i: size_t;

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

    free(f);
    return 0;
}

/**
 * Returns vertices surrounding given triangle face.
 */
ccdPtFaceVec3 :: inline (face: *ccd_pt_face_t, a: **Vector3, b: **Vector3, c: **Vector3) {
    a.* = *face.edge[0].vertex[0].v.v;
    b.* = *face.edge[0].vertex[1].v.v;

    if (face.edge[1].vertex[0] != face.edge[0].vertex[0]
            && face.edge[1].vertex[0] != face.edge[0].vertex[1]){
        c.* = *face.edge[1].vertex[0].v.v;
    }else{
        c.* = *face.edge[1].vertex[1].v.v;
    }
}

ccdPtFaceVertices :: inline (face: *ccd_pt_face_t, a: **Vector3, b: **Vector3, c: **Vector3) {
    a.* = face.edge[0].vertex[0];
    b.* = face.edge[0].vertex[1];

    if (face.edge[1].vertex[0] != face.edge[0].vertex[0]
            && face.edge[1].vertex[0] != face.edge[0].vertex[1]){
        c.* = face.edge[1].vertex[0];
    }else{
        c.* = face.edge[1].vertex[1];
    }
}

ccdPtFaceEdges :: inline (f: *ccd_pt_face_t, a: **Vector3, b: **Vector3, c: **Vector3) {
    a.* = f.edge[0];
    b.* = f.edge[1];
    c.* = f.edge[2];
}

ccdPtEdgeVec3 :: inline (e: *ccd_pt_edge_t, a: **Vector3, b: **Vector3) {
    a.* = *e.vertex[0].v.v;
    b.* = *e.vertex[1].v.v;
}

ccdPtEdgeVertices :: inline (e: *ccd_pt_edge_t, a: **Vector3, b: **Vector3) {
    a.* = e.vertex[0];
    b.* = e.vertex[1];
}

ccdPtEdgeFaces :: inline (e: *ccd_pt_edge_t, f1: **ccd_pt_face_t, f2: **ccd_pt_face_t) {
    f1.* = e.faces[0];
    f2.* = e.faces[1];
}

_ccdPtNearestUpdate :: inline (pt: *ccd_pt_t, el: *ccd_pt_el_t) {
    if (ccdEq(pt.nearest_dist, el.dist)){
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

_ccdPtNearestRenew :: (pt: *ccd_pt_t) {
    v: *ccd_pt_vertex_t;
    e: *ccd_pt_edge_t;
    f: *ccd_pt_face_t;

    pt.nearest_dist = CCD_REAL_MAX;
    pt.nearest_type = 3;
    pt.nearest = null;

    ccdListForEachEntry(*pt.vertices, v, ccd_pt_vertex_t, list, #code {
        _ccdPtNearestUpdate(pt, cast(*ccd_pt_el_t)v);
    });

    ccdListForEachEntry(*pt.edges, e, ccd_pt_edge_t, list, #code {
        _ccdPtNearestUpdate(pt, cast(*ccd_pt_el_t)e);
    });

    ccdListForEachEntry(*pt.faces, f, ccd_pt_face_t, list, #code {
        _ccdPtNearestUpdate(pt, cast(*ccd_pt_el_t)f);
    });
}



ccdPtInit :: (pt: *ccd_pt_t) {
    ccdListInit(*pt.vertices);
    ccdListInit(*pt.edges);
    ccdListInit(*pt.faces);

    pt.nearest = null;
    pt.nearest_dist = CCD_REAL_MAX;
    pt.nearest_type = 3;
}

ccdPtDestroy :: (pt: *ccd_pt_t) {
    f, f2: *ccd_pt_face_t;
    e, e2: *ccd_pt_edge_t;
    v, v2: *ccd_pt_vertex_t;

    // first delete all faces
    ccdListForEachEntrySafe(*pt.faces, f, ccd_pt_face_t, f2, ccd_pt_face_t, list, #code {
        ccdPtDelFace(pt, f);
    });

    // delete all edges
    ccdListForEachEntrySafe(*pt.edges, e, ccd_pt_edge_t, e2, ccd_pt_edge_t, list, #code {
        ccdPtDelEdge(pt, e);
    });

    // delete all vertices
    ccdListForEachEntrySafe(*pt.vertices, v, ccd_pt_vertex_t, v2, ccd_pt_vertex_t, list, #code {
        ccdPtDelVertex(pt, v);
    });
}


ccdPtAddVertex :: (pt: *ccd_pt_t, v: *ccd_support_t) -> *ccd_pt_vertex_t {
    vert: *ccd_pt_vertex_t;

    vert = CCD_ALLOC(ccd_pt_vertex_t);
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
    _ccdPtNearestUpdate(pt, cast(*ccd_pt_el_t)vert);

    return vert;
}

/**
 * Adds edge to polytope.
 */
ccdPtAddEdge :: (pt: *ccd_pt_t, v1: *ccd_pt_vertex_t, v2: *ccd_pt_vertex_t) -> *ccd_pt_edge_t {
    a, b: *Vector3;
    edge: *ccd_pt_edge_t;

    if (v1 == null || v2 == null)
        return null;

    edge = CCD_ALLOC(ccd_pt_edge_t);
    if (edge == null)
        return null;

    edge.type = CCD_PT_EDGE;
    edge.vertex[0] = v1;
    edge.vertex[1] = v2;
    edge.faces[0] = edge.faces[1] = null;

    a = *edge.vertex[0].v.v;
    b = *edge.vertex[1].v.v;
    edge.dist = ccdVec3PointSegmentDist2(ccd_vec3_origin, a, b, *edge.witness);

    ccdListAppend(*edge.vertex[0].edges, *edge.vertex_list[0]);
    ccdListAppend(*edge.vertex[1].edges, *edge.vertex_list[1]);

    ccdListAppend(*pt.edges, *edge.list);

    // update position in .nearest array
    _ccdPtNearestUpdate(pt, cast(*ccd_pt_el_t)edge);

    return edge;
}

/**
 * Adds face to polytope.
 */
ccdPtAddFace :: (pt: ccd_pt_t, e1: *ccd_pt_edge_t, e2: *ccd_pt_edge_t, e3: *ccd_pt_edge_t) -> *ccd_pt_face_t {
    a, b, c: *Vector3;
    face: *ccd_pt_face_t;
    e: *ccd_pt_edge_t;
    i: size_t;

    if (e1 == null || e2 == null || e3 == null)
        return null;

    face = CCD_ALLOC(ccd_pt_face_t);
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
    face.dist = ccdVec3PointTriDist2(ccd_vec3_origin, a, b, c, *face.witness);

    for i : 0..2 { // (i = 0; i < 3; i++){
        if (face.edge[i].faces[0] == null){
            face.edge[i].faces[0] = face;
        }else{
            face.edge[i].faces[1] = face;
        }
    }

    ccdListAppend(*pt.faces, *face.list);

    // update position in .nearest array
    _ccdPtNearestUpdate(pt, cast(*ccd_pt_el_t)face);

    return face;
}

/**
 * Recompute distances from origin for all elements in pt.
 */
ccdPtRecomputeDistances :: (pt: *ccd_pt_t) {
    v: *ccd_pt_vertex_t;
    e: *ccd_pt_edge_t;
    f: *ccd_pt_face_t;
    a, b, c: *Vector3;
    dist: ccd_real_t;

    ccdListForEachEntry(*pt.vertices, v, ccd_pt_vertex_t, list, #code {
        dist = ccdVec3Len2(*v.v.v);
        v.dist = dist;
        ccdVec3Copy(*v.witness, *v.v.v);
    });

    ccdListForEachEntry(*pt.edges, e, ccd_pt_edge_t, list, #code {
        a = *e.vertex[0].v.v;
        b = *e.vertex[1].v.v;
        dist = ccdVec3PointSegmentDist2(ccd_vec3_origin, a, b, *e.witness);
        e.dist = dist;
    });

    ccdListForEachEntry(*pt.faces, f, ccd_pt_face_t, list, #code {
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
ccdPtNearest :: (pt: *ccd_pt_t) -> *ccd_pt_el_t {
    if (!pt.nearest){
        _ccdPtNearestRenew(pt);
    }
    return pt.nearest;
}


// JAI: these are just debug functions i think so we probably dont need em rn.

// ccdPtDumpSVT :: (pt: *ccd_pt_t, fn: *char) { // fn should be what?? it was a `const char *fn`
//     FILE *fout;

//     fout = fopen(fn, "a");
//     if (fout == null)
//         return;

//     ccdPtDumpSVT2(pt, fout);

//     fclose(fout);
// }

// void ccdPtDumpSVT2(ccd_pt_t *pt, FILE *fout)
// {
//     ccd_pt_vertex_t *v, *a, *b, *c;
//     ccd_pt_edge_t *e;
//     ccd_pt_face_t *f;
//     size_t i;

//     fprintf(fout, "-----\n");

//     fprintf(fout, "Points:\n");
//     i = 0;
//     ccdListForEachEntry(&pt->vertices, v, ccd_pt_vertex_t, list){
//         v->id = i++;
//         fprintf(fout, "%lf %lf %lf\n",
//                 ccdVec3X(&v->v.v), ccdVec3Y(&v->v.v), ccdVec3Z(&v->v.v));
//     }

//     fprintf(fout, "Edges:\n");
//     ccdListForEachEntry(&pt->edges, e, ccd_pt_edge_t, list){
//         fprintf(fout, "%d %d\n", e->vertex[0]->id, e->vertex[1]->id);
//     }

//     fprintf(fout, "Faces:\n");
//     ccdListForEachEntry(&pt->faces, f, ccd_pt_face_t, list){
//         a = f->edge[0]->vertex[0];
//         b = f->edge[0]->vertex[1];
//         c = f->edge[1]->vertex[0];
//         if (c == a || c == b){
//             c = f->edge[1]->vertex[1];
//         }
//         fprintf(fout, "%d %d %d\n", a->id, b->id, c->id);
//     }
// }