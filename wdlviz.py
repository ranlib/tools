#!/usr/bin/env python3
import os
import sys
import WDL
import graphviz


def main(args):
    # load WDL document given local filename
    doc = WDL.load(args[0] if args else "/dev/stdin")
    assert doc.workflow, "No workflow in WDL document"

    # visualize workflow
    wdlviz(doc.workflow).render("workflow.dot", view=True)


def wdlviz(workflow):
    top = graphviz.Digraph(comment=workflow.name)
    top.attr(compound="true")
    node_ids = set()

    def add_node(dot, elt):
        nonlocal node_ids
        shape = None
        if isinstance(elt, WDL.WorkflowSection):
            with dot.subgraph(name="cluster-" + elt.workflow_node_id) as subdot:
                label = "scatter" if isinstance(elt, WDL.Scatter) else "if"
                subdot.attr(label=label + f"({str(elt.expr)})", rank="same")
                for child in elt.body:
                    add_node(subdot, child)
                subdot.node(
                    elt.workflow_node_id, "", style="invis", height="0", width="0", margin="0"
                )
            node_ids.add(elt.workflow_node_id)
            node_ids |= set(g.workflow_node_id for g in elt.gathers.values())
        elif isinstance(elt, WDL.Call):
            shape = "cds"
        elif isinstance(elt, WDL.Decl) and node_ids.intersection(elt.workflow_node_dependencies):
            shape = "plaintext"

        if shape:
            dot.node(elt.workflow_node_id, elt.name, shape=shape)
            node_ids.add(elt.workflow_node_id)

    for elt in workflow.body:
        add_node(top, elt)

    def add_edges(elt):
        for dep_id in elt.workflow_node_dependencies:
            dep = workflow.get_node(dep_id)
            if isinstance(dep, WDL.Tree.Gather):
                dep = dep.final_referee
                dep_id = dep.workflow_node_id
            if elt.workflow_node_id in node_ids and dep_id in node_ids:
                lhead = None
                if isinstance(elt, WDL.WorkflowSection):
                    lhead = "cluster-" + elt.workflow_node_id
                top.edge(dep_id, elt.workflow_node_id, lhead=lhead)
        if isinstance(elt, WDL.WorkflowSection):
            for child in elt.body:
                add_edges(child)

    for elt in workflow.body:
        add_edges(elt)

    return top


if __name__ == "__main__":
    main(sys.argv[1:])
