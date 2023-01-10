
def simplify_mtg(mtg, vid):
    tokeep = set([lvid for lvid in mtg.Descendants(vid) ]+[vid]+[lvid for lvid in mtg.Ancestors(vid)])
    from copy import deepcopy
    newmtg = deepcopy(mtg)
    for lvid in list(newmtg.vertices(scale=3)):
        if not lvid in tokeep and lvid in newmtg.vertices(scale=3):
            newmtg.remove_tree(lvid)
    for propvalues in newmtg.properties().values():
        for lvid in list(propvalues.keys()):
            if not lvid in tokeep:
                del propvalues[lvid]
    return newmtg
