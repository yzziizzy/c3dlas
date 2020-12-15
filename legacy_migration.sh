#!/bin/bash

# This script migrates a legacy codebase to the new naming scheme 

FN_NAMES="boxCenter|boxCenter2|boxCenter2i|boxContainsPoint|boxContainsPoint2|boxContainsPoint2i|boxDisjoint|boxDisjoint2|boxDisjoint2i|boxOverlaps|boxOverlaps2|boxOverlaps2i|boxQuadrant2|boxQuadrant2i|boxRayIntersect|boxRayIntersectFast|boxSize|boxSize2|boxSize2i|evalBezier|evalBezierNorm|evalBezierTangent|evalQBezier|evalQBezier2D|linePlaneClip|makeRay|vMatrixMul|vMatrixMulf|planeClassifyPoint|planeClassifyPointEps|planeCopy|planeFromTriangle|planeInverse|planeLineFindIntersect|planeLineFindIntersectFast|planePointDist|planePointDistSigned|pvDist|quadCenter2|quadCenterp|quadRoundInward2|quadRoundOutward2|shortestLineFromRayToRay|triArea2|triPlaneClip|triPlaneTestIntersect|triPointInside2|vAdd|vAdd2|vAdd2i|vCopy|vCopy2|vCopy2i|vCross|vDist|vDist2|vDist2i|vDistSq|vDot|vDot2|vDot2i|vEq|vEq2|vEq2i|vEqEp|vEqEp2|vInverse|vInverse2|vLerp|vLerp2|vLerp4|vMag|vMag2|vMatrixMul|vMatrixMulf|vMax|vMax2|vMax2i|vMin|vMin2|vMin2i|vNorm|vNorm2|vPointAvg|vProject|vProjectNorm|vProjectOntoPlane|vProjectOntoPlaneNormalized|vpTriFaceNormal|vRandom|vRandomNorm|vReflectAcross|vReflectAcross2|vRoundAway2|vRoundToward2|vScalarTriple|vScale|vScale2|vScale2i|vSet|vSet2|vSet2i|vSub|vSub2|vSub2i|vSwap|vSwap2|vSwap2i|vTriFaceNormal|vUnit|vUnit2"
FN3_NAMES="boxCenter|boxContainsPoint|boxDisjoint|boxOverlaps|boxRayIntersect|boxRayIntersectFast|boxSize|evalBezier|evalBezierNorm|evalBezierTangent|evalQBezier|evalQBezier2D|linePlaneClip|makeRay|vMatrixMul|vMatrixMulf|planeClassifyPoint|planeClassifyPointEps|planeCopy|planeFromTriangle|planeInverse|planeLineFindIntersect|planeLineFindIntersectFast|planePointDist|planePointDistSigned|pvDist|quadCenterp|shortestLineFromRayToRay|triPlaneClip|triPlaneTestIntersect|vAdd|vCopy|vCross|vDist|vDistSq|vDot|vEq|vEqEp|vInverse|vLerp|vMag|vMatrixMul|vMatrixMulf|vMax|vMin|vNorm|vPointAvg|vProject|vProjectNorm|vProjectOntoPlane|vProjectOntoPlaneNormalized|vpTriFaceNormal|vRandom|vRandomNorm|vReflectAcross|vScalarTriple|vScale|vSet|vSub|vSwap|vTriFaceNormal|vUnit"
ST3_NAMES="Vector|AABB|BezierSplineSegment|BezierSpline|LineSegment|Ray|Quad"

# FN_NAMES=`cat ./fns.txt | tr '\n' '|' | sed  's/\|\s*$//'` 
# echo $FN_NAMES

egrep -rlZ "\b(${FN_NAMES})" | egrep -zZ "\.[ch]" | xargs -0 \
	sed -Ei "s:\b(${FN_NAMES})\(:\1p(:g"

# cat fns.txt | egrep -v '[24]i?$' | tr '\n' '|' | sed  's/\|\s*$//'

#cat c3dlas.h |  sed -E "s:(boxSize2i|boxCenter2i):lol\1p(:g"

egrep -rlZ "\b(${FN3_NAMES}|${ST3_NAMES})" | egrep -zZ "\.[ch]" | xargs -0 \
	sed -Ei "s:\b(${FN3_NAMES})p\(:\13p(:g;s:\b(${ST3_NAMES})\b:\13:g"



