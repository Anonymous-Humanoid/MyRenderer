from myRenderer import KnotVector
from myRenderer.KnotVector import Number
import numpy as np
from typing import Tuple

def test(unnormalized_knots: Tuple[Number], normalized_knots: Tuple[Number]) -> None:
    U = KnotVector(unnormalized_knots)
    normalized_unique_knots = list(set(U.toNormalizedVector()[4:-4]))
    normalized_unique_knots.sort()
    normalized_unique_knots = tuple(normalized_unique_knots)
    
    # Testing normalization
    assert tuple(U.toNormalizedVector()[4:-4]) == normalized_knots
    assert tuple(U.toNormalizedVector()[:4]) == (0, 0, 0, 0)
    assert tuple(U.toNormalizedVector()[-4:]) == (1, 1, 1, 1)
    
    # Testing denormalization
    assert tuple(U.unnormalize(knot) for knot in U.toNormalizedVector()[4:-4]) == unnormalized_knots
    assert tuple(U.toUnnormalizedVector()[4:-4]) == unnormalized_knots
    assert tuple(U.normalize(knot) for knot in unnormalized_knots) == normalized_knots
    
    # Testing knot vector elements
    assert (U[0], U[1], U[2], U[3]) == (0, 0, 0, 0)
    assert (U[-1], U[-2], U[-3], U[-4]) == (1, 1, 1, 1)
    assert tuple(knot for knot in U) == (0,) + normalized_unique_knots + (1,)
    for i in range(4, U.total_unique_knots - 4):
        assert U[i] not in {0, 1}
    for knot1, knot2 in zip(U.toUnnormalizedVector()[:-1], U.toUnnormalizedVector()[1:]):
        assert knot2 >= knot1
    for knot1, knot2 in zip(U.toNormalizedVector()[:-1], U.toNormalizedVector()[1:]):
        assert knot2 >= knot1
    
    # Testing knot intervals
    old_k, old_i = U.first_outer_knot_mul, 1
    assert U.getKnotIntervalIndices(-np.inf) is None
    assert U.getKnotIntervalIndices(np.inf) == (len(U) - 1, U.total_unique_knots - 1)
    for knot in U.getUniqueKnots():
        indices = U.getKnotIntervalIndices(knot)
        if knot <= KnotVector.FIRST_NORMALIZED_OUTER_KNOT:
            assert indices is None
        else:
            k, i = indices
            print(k, i, '-', old_k, old_i)
            if k >= len(U) and i >= U.total_unique_knots:
                assert indices == (len(U) - 1, U.total_unique_knots - 1)
            # TODO Fix assertions
            elif k > old_k and i > old_i:
                assert i == old_i + 1
                assert k == old_k + U.multiplicity(U.unnormalize(U[k]))
            else:
                assert indices == (old_k, old_i)
            old_k, old_i = indices

if __name__ == '__main__':
    test((5, 6, 7, 8), (1/5, 2/5, 3/5, 4/5))
    test((1, 1, 1, 2, 2, 2), (1/3, 1/3, 1/3, 2/3, 2/3, 2/3))
    test((1, 2, 3, 4, 5, 6, 7), (1/8, 2/8, 3/8, 4/8, 5/8, 6/8, 7/8))
    test((1, 7), (1/8, 7/8))
    test((0, 1), (1/3, 2/3))
    print('All tests passed!')