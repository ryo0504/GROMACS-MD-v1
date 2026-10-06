def rdf(universe, selection_a, selection_b, bins=75, distance_range=(0.0, 15.0)):
    """Distances in Angstrom; explicit selections must define self/exclusion policy."""
    from MDAnalysis.analysis.rdf import InterRDF
    result = InterRDF(universe.select_atoms(selection_a), universe.select_atoms(selection_b),
                      nbins=bins, range=distance_range).run()
    return result.results.bins, result.results.rdf
