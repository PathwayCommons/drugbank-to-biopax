package com.google.gsoc14.drugbank2biopax.controller;

import ca.drugbank.model.*;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.biopax.paxtools.model.BioPAXElement;
import org.biopax.paxtools.model.BioPAXFactory;
import org.biopax.paxtools.model.BioPAXLevel;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.model.level3.*;
import org.biopax.paxtools.model.level3.Process;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBElement;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Unmarshaller;
import java.io.InputStream;
import java.util.List;

public class DrugBank2BioPAXConverter {
    private static Log log = LogFactory.getLog(DrugBank2BioPAXConverter.class);

    private BioPAXFactory bioPAXFactory = BioPAXLevel.L3.getDefaultFactory();

    public BioPAXFactory getBioPAXFactory() {
        return bioPAXFactory;
    }

    public void setBioPAXFactory(BioPAXFactory bioPAXFactory) {
        this.bioPAXFactory = bioPAXFactory;
    }

    protected <T extends BioPAXElement> T create(Class<T> aClass, String id) {
        return getBioPAXFactory().create(aClass, id);
    }

    public Model convert(InputStream drugBankStream) throws JAXBException {
        Model model = BioPAXLevel.L3.getDefaultFactory().createModel();

        JAXBContext jaxbContext = JAXBContext.newInstance("ca.drugbank.model");
        Unmarshaller unmarshaller = jaxbContext.createUnmarshaller();
        Object unmarshalled = unmarshaller.unmarshal(drugBankStream);
        DrugbankType value = (DrugbankType) ((JAXBElement) unmarshalled).getValue();
        List<DrugType> drugs = value.getDrug();
        log.info(drugs.size() + " drugs in the database.");

        int drugTargets = 0;
        for (DrugType drug : drugs) {
            SmallMolecule smallMolecule = convertDrugToSmallMolecule(model, drug);
            for (TargetType targetType : drug.getTargets().getTarget()) {
                convertTargetToReaction(model, smallMolecule, targetType);
                drugTargets++;
            }
        }

        log.info("Conversion done. " + drugTargets + " drug targets converted into reactions.");

        return model;
    }

    private Control convertTargetToReaction(Model model, SmallMolecule smallMolecule, TargetType targetType) {
        Control control = create(Control.class, "control_" + targetType.hashCode());
        model.add(control);
        control.addController(smallMolecule);
        for (PolypeptideType polypeptideType : targetType.getPolypeptide()) {
            control.addControlled(createInhibitionReaction(model, polypeptideType));
        }
        control.setControlType(decideOnControlType(targetType));
        String name = constructName(smallMolecule, targetType);
        if(name != null) {
            control.setStandardName(name);
            control.setDisplayName(name);
            control.addName(name);
        }

        return control;
    }

    private ControlType decideOnControlType(TargetType targetType) {
        // this is the reaction that inhibits the target. So the control type should be the reverse of the term
        for (String axn : targetType.getActions().getAction()) {
            if(axn.equalsIgnoreCase("substrate") ||
                    axn.equalsIgnoreCase("agonist") ||
                    axn.equalsIgnoreCase("inducer") ||
                    axn.equalsIgnoreCase("potentiator") ||
                    axn.equalsIgnoreCase("stimulator") ||
                    axn.equalsIgnoreCase("cofactor") ||
                    axn.equalsIgnoreCase("ligand") )
            {
                return ControlType.INHIBITION;
            }
        }
        return ControlType.ACTIVATION;
    }

    private String constructName(SmallMolecule smallMolecule, TargetType targetType) {
        if(targetType.getKnownAction().value().equalsIgnoreCase("yes")) {
            String actionStr = "";
            for (String action : targetType.getActions().getAction()) {
                actionStr += action + ", ";
            }

            if(actionStr.isEmpty()) {
                return null;
            } else {
                actionStr = actionStr.substring(0, actionStr.length() - 2);
                return smallMolecule.getDisplayName() + " acts as " + actionStr;
            }
        }

        return null;
    }

    private Process createInhibitionReaction(Model model, PolypeptideType polypeptideType) {
        String pid = polypeptideType.getId();
        String rxnId = "rxn_" + pid;
        BiochemicalReaction rxn = (BiochemicalReaction) model.getByID(rxnId);
        if(rxn != null) {
            return rxn;
        }

        rxn = create(BiochemicalReaction.class, rxnId);
        model.add(rxn);
        rxn.addLeft(createSPEFromPolypeptide(model, polypeptideType, "active"));
        rxn.addRight(createSPEFromPolypeptide(model, polypeptideType, "inactive"));
        rxn.setConversionDirection(ConversionDirectionType.LEFT_TO_RIGHT);

        return rxn;
    }

    private PhysicalEntity createSPEFromPolypeptide(Model model, PolypeptideType polypeptideType, String term) {
        String pid = polypeptideType.getId();
        String refId = "ref_" + pid;
        String speId = term + "_" + pid;
        ProteinReference proteinReference = (ProteinReference) model.getByID(refId);
        if(proteinReference == null) {
            proteinReference = create(ProteinReference.class, refId);
            model.add(proteinReference);
            setPolypeptideNames(polypeptideType, proteinReference);

            UnificationXref unificationXref = create(UnificationXref.class, "uxref_" + pid);
            model.add(unificationXref);
            String source = polypeptideType.getSource();
            if(source == null || source.isEmpty()) {
                source = "UniProt KnowledgeBase";
            }
            unificationXref.setDb(source);
            unificationXref.setId(pid);
            proteinReference.addXref(unificationXref);

            for (PolypeptideExternalIdentifierType externalIdentifierType : polypeptideType.getExternalIdentifiers().getExternalIdentifier()) {
                String exid = externalIdentifierType.getIdentifier();
                String exdb = externalIdentifierType.getResource().value();
                RelationshipXref relationshipXref = create(RelationshipXref.class, "rxref_" + externalIdentifierType.hashCode());
                model.add(relationshipXref);
                relationshipXref.setDb(exdb);
                relationshipXref.setId(exid);
                proteinReference.addXref(relationshipXref);
            }

            String organism = polypeptideType.getOrganism().getValue();
            if(organism != null && !organism.isEmpty()) {
                String bioSrcId = "biosrc_" + organism.hashCode();
                BioSource bioSource = (BioSource) model.getByID(bioSrcId);
                if (bioSource == null) {
                    bioSource = create(BioSource.class, bioSrcId);
                    model.add(bioSource);
                    bioSource.setStandardName(organism);
                    bioSource.setDisplayName(organism);
                    bioSource.addName(organism);
                }
                proteinReference.setOrganism(bioSource);
            }
        }

        Protein protein = (Protein) model.getByID(speId);
        if(protein == null) {
            protein = create(Protein.class, speId);
            model.add(protein);
            setPolypeptideNames(polypeptideType, protein);
            protein.setEntityReference(proteinReference);

            ModificationFeature modificationFeature = create(ModificationFeature.class, "mod_" + term + "_" + pid);
            model.add(modificationFeature);

            modificationFeature.setModificationType(createSeqModVocab(model, term));
            protein.addFeature(modificationFeature);
            protein.getEntityReference().addEntityFeature(modificationFeature);
        }

        return protein;
    }

    private SequenceModificationVocabulary createSeqModVocab(Model model, String term) {
        String vocabId = "vocab_" + term;
        SequenceModificationVocabulary sequenceModificationVocabulary = (SequenceModificationVocabulary) model.getByID(vocabId);
        if(sequenceModificationVocabulary == null) {
            sequenceModificationVocabulary = create(SequenceModificationVocabulary.class, vocabId);
            model.add(sequenceModificationVocabulary);
            sequenceModificationVocabulary.addTerm(term);
        }

        return sequenceModificationVocabulary;
    }

    private SmallMolecule convertDrugToSmallMolecule(Model model, DrugType drug) {
        String rdfId = createRDFId(drug);
        SmallMoleculeReference reference = create(SmallMoleculeReference.class, "ref_" + rdfId);
        model.add(reference);
        SmallMolecule smallMolecule = create(SmallMolecule.class, rdfId);
        model.add(smallMolecule);
        smallMolecule.setEntityReference(reference);

        setDrugNames(drug, smallMolecule);
        setDrugNames(drug, reference);

        String primaryDrugId = getPrimaryDrugId(drug);
        UnificationXref unificationXref = create(UnificationXref.class, "uxref_" + primaryDrugId);
        model.add(unificationXref);
        unificationXref.setDb("DrugBank");
        unificationXref.setId(primaryDrugId);
        reference.addXref(unificationXref);

        for (ExternalIdentifierType externalIdentifierType : drug.getExternalIdentifiers().getExternalIdentifier()) {
            RelationshipXref relationshipXref = create(RelationshipXref.class, "rxref_" + externalIdentifierType.hashCode());
            model.add(relationshipXref);
            relationshipXref.setId(externalIdentifierType.getIdentifier());
            relationshipXref.setDb(externalIdentifierType.getResource().value());
            reference.addXref(relationshipXref);
        }

        smallMolecule.addComment(drug.getDescription());
        smallMolecule.addComment(drug.getGeneralReferences());

        for (ExperimentalPropertyType experimentalPropertyType : drug.getExperimentalProperties().getProperty()) {
            String kind = experimentalPropertyType.getKind().value();
            String value = experimentalPropertyType.getValue();
            if(kind.equalsIgnoreCase("Molecular Weight")) {
                try {
                    reference.setMolecularWeight(Float.parseFloat(value));
                } catch(NumberFormatException e) {
                    log.warn("Couldn't parse molecular weight for drug " + primaryDrugId + ": " + value);
                }
            } else if(kind.equalsIgnoreCase("Molecular Formula")) {
                reference.setChemicalFormula(value);
            } else if(kind.equalsIgnoreCase("SMILES")) {
                ChemicalStructure chemicalStructure = create(ChemicalStructure.class, "smiles_" + experimentalPropertyType.hashCode());
                model.add(chemicalStructure);
                chemicalStructure.setStructureFormat(StructureFormatType.SMILES);
                chemicalStructure.setStructureData(value);
            }
        }

        // TODO: create sub reactions for this drug (chemical -> chemical)

        return smallMolecule;
    }

    private void setPolypeptideNames(PolypeptideType polypeptideType, Named named) {
        String name = polypeptideType.getGeneName();
        if(name == null || name.isEmpty()) {
            name = polypeptideType.getName();
        } else {
            named.addName(polypeptideType.getName());
        }

        named.setDisplayName(name);
        named.setStandardName(name);
        named.addName(name);

        for (String s : polypeptideType.getSynonyms().getSynonym()) {
            named.addName(s);
        }
    }


    private void setDrugNames(DrugType drug, Named named) {
        String name = drug.getName();

        named.setDisplayName(name);
        named.setStandardName(name);
        named.addName(name);

        for (SynonymType synonymType : drug.getSynonyms().getSynonym()) {
            named.addName(synonymType.getValue());
        }

        for (BrandType brandType : drug.getBrands().getBrand()) {
            named.addName(brandType.getValue());
        }
    }

    private String createRDFId(DrugType drug) {
        String primaryDrugId = getPrimaryDrugId(drug);
        return primaryDrugId != null ? primaryDrugId : drug.hashCode() + "";
    }

    private String getPrimaryDrugId(DrugType drug) {
        for (DrugbankDrugIdType drugbankDrugIdType : drug.getDrugbankId()) {
            if(drugbankDrugIdType.isPrimary()) {
                return drugbankDrugIdType.getValue();
            }
        }

        return null;
    }

}
