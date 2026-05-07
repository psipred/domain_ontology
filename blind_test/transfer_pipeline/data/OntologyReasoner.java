package org.ontology;

import org.semanticweb.HermiT.ReasonerFactory;
import org.semanticweb.owlapi.apibinding.OWLManager;
import org.semanticweb.owlapi.model.*;
import org.semanticweb.owlapi.reasoner.*;
import org.semanticweb.owlapi.util.InferredSubClassAxiomGenerator;
import org.semanticweb.owlapi.util.InferredOntologyGenerator;

import java.io.*;
import java.util.*;
import java.util.stream.Collectors;

/**
 * OntologyReasoner.java
 * ─────────────────────────────────────────────────────────────────────────────
 * STAGE 4 — OWL Reasoning with OWLAPI + HermiT
 *
 * Loads core_ontology.owl, checks consistency, lists unsatisfiable classes,
 * and prints inferred subclass hierarchies for key domain ontology classes.
 *
 * Usage:
 *   java -jar ontology-reasoner.jar <ontology_path> [output_report_path]
 *
 * Build:
 *   cd java_reasoner && mvn clean package
 *   java -jar target/ontology-reasoner-1.0.jar input/core_ontology.owl output/reasoner_report.txt
 *
 * ADAPT: Update TARGET_CLASSES to match your ontology class IRIs.
 * ─────────────────────────────────────────────────────────────────────────────
 */
public class OntologyReasoner {

    // ── ADAPT: update these to your ontology namespace ────────────────────────
    private static final String ONTOLOGY_IRI_BASE =
        "http://example.org/domain-ontology-core#";

    // Classes to inspect for inferred subclass hierarchies
    private static final List<String> TARGET_CLASSES = Arrays.asList(
        "Protein",
        "DomainType",
        "DomainOccurrence",
        "PositionalCategoryTerm",
        "FunctionTagTerm",
        "BindingTargetTerm",
        "GOTerm",
        "ReactomePathway",
        "KEGGEntry",
        "CATHSuperfamily"
    );

    public static void main(String[] args) throws Exception {
        String ontologyPath  = args.length > 0 ? args[0] : "input/core_ontology.owl";
        String reportPath    = args.length > 1 ? args[1] : "output/reasoner_report.txt";

        System.out.println("======================================================");
        System.out.println("  cAMP Pathway Domain Ontology — OWL Reasoner");
        System.out.println("======================================================");
        System.out.println("  Ontology: " + ontologyPath);
        System.out.println("  Report:   " + reportPath);

        // ── Load ontology ─────────────────────────────────────────────────────
        OWLOntologyManager manager = OWLManager.createOWLOntologyManager();
        File ontFile = new File(ontologyPath);
        if (!ontFile.exists()) {
            System.err.println("[ERROR] Ontology file not found: " + ontologyPath);
            System.exit(1);
        }
        OWLOntology ontology = manager.loadOntologyFromOntologyDocument(ontFile);
        OWLDataFactory factory = manager.getOWLDataFactory();
        System.out.println("  Loaded " + ontology.getAxiomCount() + " axioms.");

        // ── Create HermiT reasoner ────────────────────────────────────────────
        System.out.println("\nStarting HermiT reasoner ...");
        OWLReasonerFactory reasonerFactory = new ReasonerFactory();
        ConsoleProgressMonitor progressMonitor = new ConsoleProgressMonitor();
        OWLReasonerConfiguration config = new SimpleConfiguration(progressMonitor);
        OWLReasoner reasoner = reasonerFactory.createReasoner(ontology, config);
        reasoner.precomputeInferences(
            InferenceType.CLASS_HIERARCHY,
            InferenceType.CLASS_ASSERTIONS,
            InferenceType.OBJECT_PROPERTY_HIERARCHY
        );

        PrintWriter out = new PrintWriter(new FileWriter(reportPath));
        StringBuilder sb = new StringBuilder();

        // ── Consistency check ─────────────────────────────────────────────────
        boolean consistent = reasoner.isConsistent();
        sb.append("=== ONTOLOGY CONSISTENCY ===\n");
        sb.append("Consistent: ").append(consistent).append("\n\n");
        System.out.println("  Consistent: " + consistent);

        if (!consistent) {
            sb.append("[ERROR] Ontology is INCONSISTENT.\n");
            sb.append("Cannot compute further inferences.\n");
            out.print(sb);
            out.close();
            reasoner.dispose();
            return;
        }

        // ── Unsatisfiable classes ─────────────────────────────────────────────
        sb.append("=== UNSATISFIABLE CLASSES ===\n");
        Node<OWLClass> bottomNode = reasoner.getUnsatisfiableClasses();
        Set<OWLClass> unsatisfiable = bottomNode.getEntitiesMinusBottom();
        if (unsatisfiable.isEmpty()) {
            sb.append("None (all classes are satisfiable).\n");
        } else {
            for (OWLClass cls : unsatisfiable) {
                sb.append("  UNSATISFIABLE: ").append(getLocalName(cls)).append("\n");
                System.out.println("  [WARN] Unsatisfiable: " + getLocalName(cls));
            }
        }
        sb.append("\n");

        // ── Ontology statistics ───────────────────────────────────────────────
        sb.append("=== ONTOLOGY STATISTICS ===\n");
        long nClasses = ontology.classesInSignature().count();
        long nObjProp = ontology.objectPropertiesInSignature().count();
        long nDataProp = ontology.dataPropertiesInSignature().count();
        long nIndiv   = ontology.individualsInSignature().count();
        sb.append("Classes:           ").append(nClasses).append("\n");
        sb.append("Object properties: ").append(nObjProp).append("\n");
        sb.append("Data properties:   ").append(nDataProp).append("\n");
        sb.append("Individuals:       ").append(nIndiv).append("\n\n");

        // ── Inferred subclass hierarchies ─────────────────────────────────────
        sb.append("=== INFERRED SUBCLASS HIERARCHIES ===\n\n");
        for (String className : TARGET_CLASSES) {
            IRI classIRI = IRI.create(ONTOLOGY_IRI_BASE + className);
            OWLClass targetClass = factory.getOWLClass(classIRI);
            if (!ontology.containsClassInSignature(classIRI)) {
                sb.append("[").append(className).append("] — NOT FOUND in ontology\n\n");
                continue;
            }
            sb.append("[").append(className).append("]\n");

            // Direct subclasses
            NodeSet<OWLClass> directSubs = reasoner.getSubClasses(targetClass, true);
            List<String> directNames = directSubs.entities()
                .filter(c -> !c.isBottomEntity())
                .map(OntologyReasoner::getLocalName)
                .sorted()
                .collect(Collectors.toList());
            sb.append("  Direct subclasses (").append(directNames.size()).append("):\n");
            for (String sub : directNames) {
                sb.append("    - ").append(sub).append("\n");
            }

            // All descendant subclasses
            NodeSet<OWLClass> allSubs = reasoner.getSubClasses(targetClass, false);
            long allCount = allSubs.entities()
                .filter(c -> !c.isBottomEntity())
                .count();
            sb.append("  All descendants: ").append(allCount).append("\n\n");
        }

        // ── Inferred object property domains/ranges ────────────────────────
        sb.append("=== OBJECT PROPERTY ASSERTIONS (sample) ===\n");
        ontology.objectPropertiesInSignature()
            .limit(20)
            .forEach(prop -> {
                Set<OWLClassExpression> domains = ontology
                    .getObjectPropertyDomainAxioms(prop).stream()
                    .map(OWLObjectPropertyDomainAxiom::getDomain)
                    .collect(java.util.stream.Collectors.toSet());
                Set<OWLClassExpression> ranges = ontology
                    .getObjectPropertyRangeAxioms(prop).stream()
                    .map(OWLObjectPropertyRangeAxiom::getRange)
                    .collect(java.util.stream.Collectors.toSet());
                String propName = getLocalName(prop);
                if (!propName.isEmpty()) {
                    sb.append("  ").append(propName).append("\n");
                    domains.forEach(d -> sb.append("    domain: ").append(d).append("\n"));
                    ranges.forEach(r  -> sb.append("    range: ").append(r).append("\n"));
                }
            });
        sb.append("\n");

        // ── Write report ──────────────────────────────────────────────────────
        out.print(sb);
        out.close();
        reasoner.dispose();

        System.out.println("\n[OK] Reasoner report written to: " + reportPath);
        System.out.println("  Unsatisfiable classes: " + unsatisfiable.size());
    }

    /** Extract the local name (fragment) from an OWL entity's IRI. */
    private static String getLocalName(OWLEntity entity) {
        String iri = entity.getIRI().toString();
        int idx = Math.max(iri.lastIndexOf('#'), iri.lastIndexOf('/'));
        return idx >= 0 ? iri.substring(idx + 1) : iri;
    }
}
