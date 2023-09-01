use super::shape;

// https://web.archive.org/web/20130126163405/http://geomalgorithms.com/a03-_inclusion.html
// Copyright 2000 softSurfer, 2012 Dan Sunday
// This code may be freely used and modified for any purpose
// providing that this copyright notice is included with it.
// SoftSurfer makes no warranty for this code, and cannot be held
// liable for any real or imagined damage resulting from its use.
// Users of this code must verify correctness for their application.

// a Point is defined by its coordinates {int x, y;}
//===================================================================

// isLeft(): tests if a point is Left|On|Right of an infinite line.
//    Input:  three points P0, P1, and P2
//    Return: >0 for P2 left of the line through P0 and P1
//            =0 for P2  on the line
//            <0 for P2  right of the line
//    See: Algorithm 1 "Area of Triangles and Polygons"
fn is_left(p0: &[f32; 2], p1: &[f32; 2], p2: &[f32; 2]) -> f32 {
    return (p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1]);
}
//===================================================================

#[inline(always)]
fn cn_pn_poly_side(v0: &[f32; 2], v1: &[f32; 2], p: &[f32; 2]) -> i32 {
    if ((v0[1] <= p[1]) && (v1[1] > p[1]))     // an upward crossing
      || ((v0[1] > p[1]) && (v1[1] <= p[1]))
    {
        // a downward crossing
        // compute the actual edge-ray intersect x-coordinate
        let vt = (p[1] - v0[1]) / (v1[1] - v0[1]);
        if p[0] < v0[0] + vt * (v1[0] - v0[0]) {
            // P[0] < intersect
            return 1; // a valid crossing of y=P[1] right of P[0]
        }
    }
    return 0;
}

// cn_PnPoly(): crossing number test for a point in a polygon
//      Input:   P = a point,
//               V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
//      Return:  0 = outside, 1 = inside
// This code is patterned after [Franklin, 2000]
pub fn cn_pn_poly(p: &[f32; 2], shape: &shape::Shape) -> i32 {
    let mut cn: i32 = 0; // the crossing number counter
                         // loop through all edges of the polygon
    for i in 0..shape.pts.len() - 1 {
        cn += cn_pn_poly_side(&shape.pts[i], &shape.pts[i + 1], &p);
    }
    cn += cn_pn_poly_side(&shape.pts[shape.pts.len() - 1], &shape.pts[0], &p);
    return cn % 2;
}

//===================================================================
#[inline(always)]
fn wn_pn_poly_side(v0: &[f32; 2], v1: &[f32; 2], p: &[f32; 2]) -> i32 {
    if v0[1] <= p[1] {
        // start y <= p[1]
        if v1[1] > p[1] {
            // an upward crossing
            if is_left(v0, v1, p) > 0. {
                // p left of  edge
                return 1; // have  a valid up intersect
            }
        }
    } else {
        // start y > p[1] (no test needed)
        if v1[1] <= p[1] {
            // a downward crossing
            if is_left(v0, v1, p) < 0. {
                // p right of  edge
                return -1; // have  a valid down intersect
            }
        }
    }
    return 0;
}

// wn_PnPoly(): winding number test for a point in a polygon
//      Input:   P = a point,
//               V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
//      Return:  wn = the winding number (=0 only when P is outside)
pub fn wn_pn_poly(p: &[f32; 2], shape: &shape::Shape) -> i32 {
    let mut wn: i32 = 0; // the  winding number counter

    for i in 0..shape.pts.len() - 1 {
        wn += wn_pn_poly_side(&shape.pts[i], &shape.pts[i + 1], &p);
    }
    wn += wn_pn_poly_side(&shape.pts[shape.pts.len() - 1], &shape.pts[0], &p);
    return wn;
}
//===================================================================