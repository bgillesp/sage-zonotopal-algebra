class SetUtils(object):
    @staticmethod
    def minimal_sets(subsets):
        """
        Returns the minimal elements (with respect to inclusion) of a collection of subsets.
        """
        def geq(s1, s2):
            return set(s1) >= set(s2)
        # pass geq to retrieve minimal sets
        return SetUtils.max_elements(subsets, geq)

    @staticmethod
    def max_elements(elts, leq):
        """
        Returns the maximal elements with respect to a less than or equal to function, leq.
        Can return minimal elements by instead passing a greater than or equal to function.
        """
        max_elts = []
        for e1 in elts:
            maximal = True
            non_max = []
            for e2 in max_elts:
                if leq(e1, e2):
                    maximal = False
                elif leq(e2, e1):
                    non_max.append(e2)
            for e2 in non_max:
                max_elts.remove(e2)
            if maximal:
                max_elts.append(e1)
        return max_elts
