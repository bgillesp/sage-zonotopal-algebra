class SetUtils(object):
    @staticmethod
    def minimal_sets(subsets):
        """
        Returns the minimal elements (with respect to inclusion) of a collection of subsets.
        """
        min_sets = []
        for s1 in subsets:
            minimal = True
            larger = []
            for s2 in min_sets:
                if set(s1) >= set(s2):
                    minimal = False
                elif set(s1) < set(s2):
                    larger.append(s2)
            for s in larger:
                min_sets.remove(s)
            if minimal:
                min_sets.append(s1)
        return min_sets
