package ca.drugbank.converter;

import ca.drugbank.model.*;
import org.apache.commons.lang3.StringUtils;
import org.biopax.paxtools.controller.ModelUtils;
import org.biopax.paxtools.model.BioPAXElement;
import org.biopax.paxtools.model.BioPAXFactory;
import org.biopax.paxtools.model.BioPAXLevel;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.model.level3.*;
import org.biopax.paxtools.model.level3.Process;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBElement;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Unmarshaller;
import java.io.InputStream;
import java.util.List;
import java.util.UUID;

public class DrugbankToBiopaxConverter {
    private final static Logger log = LoggerFactory.getLogger(DrugbankToBiopaxConverter.class);

    private final static BioPAXFactory bioPAXFactory = BioPAXLevel.L3.getDefaultFactory();

    private <T extends BioPAXElement> T create(Class<T> aClass, String id) {
        return bioPAXFactory.create(aClass, completeId(id));
    }

    private <T extends Xref> T findOrCreateXref(Model model, Class<T> xrefType, String id, String db)
    {
        //some quick-fix, normalization
        if("uniprot accession".equalsIgnoreCase(db))
            db = "UniProt"; //it's UniProt ID (accession is like Q12345, P43210)!

        String localUri = xrefType.getSimpleName() + "_" + ModelUtils.md5hex(db + id);
        T x = (T) model.getByID(completeId(localUri));
        if(x == null) {
            x = create(xrefType, localUri);
            x.setDb(db);
            x.setId(id);
            model.add(x);
        }
        return x;
    }

    private String xmlBase = "";

    public String getXmlBase() {
        return xmlBase;
    }

    public void setXmlBase(String xmlBase) {
        this.xmlBase = xmlBase;
    }

    private String completeId(String id) {
        return (StringUtils.containsAny(id,"identifiers.org/", "bioregistry.io/"))
            ? id : getXmlBase() + id;
    }

    public Model convert(InputStream drugBankStream) throws JAXBException {
        Model model = bioPAXFactory.createModel();
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

        model.setXmlBase(getXmlBase());
        return model;
    }

    private Control convertTargetToReaction(Model model, SmallMolecule smallMolecule, TargetType targetType) {
        Control control = create(Control.class, "control_" + UUID.randomUUID());
        model.add(control);
        control.addController(smallMolecule);
        for (PolypeptideType polypeptideType : targetType.getPolypeptide()) {
            if(isAcetylation(targetType)) {
                control.addControlled(createAcetylationReaction(model, polypeptideType));
            } else {
                control.addControlled(createInhibitionReaction(model, polypeptideType));
            }
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

    private boolean isAcetylation(TargetType targetType) {
        for (String axn : targetType.getActions().getAction()) {
            if(axn.equalsIgnoreCase("acetylation")) {
                return true;
            }
        }

        return false;
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
                    axn.equalsIgnoreCase("activator") ||
                    axn.equalsIgnoreCase("ligand") )
            {
                return ControlType.INHIBITION;
            }
        }
        return ControlType.ACTIVATION;
    }

    private String constructName(SmallMolecule smallMolecule, TargetType targetType) {
        if(targetType.getKnownAction().value().equalsIgnoreCase("yes")) {
            StringBuilder sb = new StringBuilder();
            for (String action : targetType.getActions().getAction()) {
                sb.append(action).append(", ");
            }
            String act = sb.toString();
            if(act.isEmpty()) {
                return null;
            } else {
                act = act.substring(0, act.length() - 2);
                return smallMolecule.getDisplayName() + " acts as " + act;
            }
        }

        return null;
    }

    private Process createInhibitionReaction(Model model, PolypeptideType polypeptideType) {
        String pid = polypeptideType.getId();
        String rxnId = "rxn_" + pid;
        BiochemicalReaction rxn = (BiochemicalReaction) model.getByID(completeId(rxnId));
        if(rxn == null) {
            rxn = create(BiochemicalReaction.class, rxnId);
            model.add(rxn);
            rxn.addLeft(createSPEFromPolypeptide(model, polypeptideType, "active"));
            rxn.addRight(createSPEFromPolypeptide(model, polypeptideType, "inactive"));
            rxn.setConversionDirection(ConversionDirectionType.LEFT_TO_RIGHT);
        }
        return rxn;
    }

    private Process createAcetylationReaction(Model model, PolypeptideType polypeptideType) {
        String pid = polypeptideType.getId();
        String rxnId = "rxn_ace_" + pid;
        BiochemicalReaction rxn = (BiochemicalReaction) model.getByID(completeId(rxnId));
        if(rxn == null) {
            rxn = create(BiochemicalReaction.class, rxnId);
            model.add(rxn);
            rxn.addLeft(createSPEFromPolypeptide(model, polypeptideType, null));
            rxn.addRight(createSPEFromPolypeptide(model, polypeptideType, "acetylated"));
            rxn.setConversionDirection(ConversionDirectionType.LEFT_TO_RIGHT);
        }
        return rxn;
    }

    private PhysicalEntity createSPEFromPolypeptide(Model model, PolypeptideType polypeptideType, String term) {
        String pid = polypeptideType.getId();
        String refId = "ref_" + pid;
        String speId = (term == null ? "nomod" : term)  + "_" + pid;
        ProteinReference proteinReference = (ProteinReference) model.getByID(completeId(refId));
        if(proteinReference == null) {
            proteinReference = create(ProteinReference.class, refId);
            model.add(proteinReference);
            setPolypeptideNames(polypeptideType, proteinReference);

            String source = polypeptideType.getSource();
            if(source == null || source.isEmpty()) {
                source = "uniprot";
            }
            UnificationXref unificationXref = findOrCreateXref(model, UnificationXref.class, pid, source);
            proteinReference.addXref(unificationXref);

            for (PolypeptideExternalIdentifierType externalIdentifierType : polypeptideType.getExternalIdentifiers().getExternalIdentifier()) {
                String exid = externalIdentifierType.getIdentifier();
                String exdb = externalIdentifierType.getResource().value();
                RelationshipXref relationshipXref = findOrCreateXref(model, RelationshipXref.class, exid, exdb);
                proteinReference.addXref(relationshipXref);
            }

            String organism = polypeptideType.getOrganism().getValue();
            if(organism != null && !organism.isEmpty()) {
                String bioSrcId = "biosrc_" + ModelUtils.md5hex(organism);
                BioSource bioSource = (BioSource) model.getByID(completeId(bioSrcId));
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

        Protein protein = (Protein) model.getByID(completeId(speId));
        if(protein == null) {
            protein = create(Protein.class, speId);
            model.add(protein);
            setPolypeptideNames(polypeptideType, protein);
            protein.setEntityReference(proteinReference);

            if(term != null) {
                ModificationFeature modificationFeature = create(ModificationFeature.class, "mod_" + term + "_" + pid);
                model.add(modificationFeature);
                modificationFeature.setModificationType(createSeqModVocab(model, term));
                protein.addFeature(modificationFeature);
                protein.getEntityReference().addEntityFeature(modificationFeature);
            }
        }

        return protein;
    }

    private SequenceModificationVocabulary createSeqModVocab(Model model, String term) {
        String vocabId = "vocab_" + term;
        SequenceModificationVocabulary sequenceModificationVocabulary = (SequenceModificationVocabulary) model.getByID(completeId(vocabId));
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
        UnificationXref unificationXref = findOrCreateXref(model, UnificationXref.class, primaryDrugId, "DrugBank");
        reference.addXref(unificationXref);

        // find/generate relationship xrefs
        for (ExternalIdentifierType externalIdentifierType : drug.getExternalIdentifiers().getExternalIdentifier()) {
            reference.addXref(findOrCreateXref(model, RelationshipXref.class,
              externalIdentifierType.getIdentifier(), externalIdentifierType.getResource().value()));
        }

        smallMolecule.addComment(drug.getDescription());
        if(drug.getGeneralReferences() != null) {
            if (drug.getGeneralReferences().getArticles() != null) {
                for (ArticleType article : drug.getGeneralReferences().getArticles().getArticle()) {
                    String uri = "bioregistry.io/pubmed:" + article.getPubmedId();
                    PublicationXref px = (PublicationXref) model.getByID(uri);
                    if(px == null) {
                        px = create(PublicationXref.class, uri);
                        model.add(px);
                        px.setId(article.getPubmedId());
                        px.setDb("pubmed");
                        px.addUrl(px.getUri());
                        px.addComment(article.getCitation());
                    }
                    reference.addXref(px);
                }
            }
            if (drug.getGeneralReferences().getLinks() != null) {
                for (LinkType link : drug.getGeneralReferences().getLinks().getLink()) {
                    String uri = "xlink_" + ModelUtils.md5hex(link.getUrl());
                    PublicationXref px = (PublicationXref) model.getByID(completeId(uri));
                    if(px == null) {
                        px = create(PublicationXref.class, uri);
                        px.addUrl(link.getUrl());
                        if(!link.getTitle().equalsIgnoreCase("link")) {
                            px.addSource(link.getTitle());
                        }
                        model.add(px);
                    }
                    reference.addXref(px);
                }
            }
            if (drug.getGeneralReferences().getTextbooks() != null) {
                for (TextbookType tb : drug.getGeneralReferences().getTextbooks().getTextbook()) {
                    String uri = "bioregistry.io/isbn:" + tb.getIsbn();
                    PublicationXref px = (PublicationXref) model.getByID(uri);
                    if(px == null) {
                        px = create(PublicationXref.class, uri);
                        model.add(px);
                        px.addUrl(uri);
                        px.setId(tb.getIsbn());
                        px.setDb("isbn");
                        px.addComment(tb.getCitation());
                    }
                    reference.addXref(px);
                }
            }
        }

        for (ExperimentalPropertyType experimentalPropertyType : drug.getExperimentalProperties().getProperty()) {
            if(experimentalPropertyType.getKind() == null)
                continue;
            String kind = experimentalPropertyType.getKind().value();
            String value = experimentalPropertyType.getValue();
            if("Molecular Weight".equalsIgnoreCase(kind)) {
                try {
                    reference.setMolecularWeight(Float.parseFloat(value));
                } catch(NumberFormatException e) {
                    log.warn("Couldn't parse molecular weight for drug " + primaryDrugId + ": " + value);
                }
            } else if(kind.equalsIgnoreCase("Molecular Formula")) {
                reference.setChemicalFormula(value);
            } else if(kind.equalsIgnoreCase("SMILES")) {
                String localId = "smiles_" + ModelUtils.md5hex(value);
                ChemicalStructure chemicalStructure = (ChemicalStructure) model.getByID(completeId(localId));
                if(chemicalStructure==null) {
                    chemicalStructure = create(ChemicalStructure.class, localId);
                    model.add(chemicalStructure);
                    chemicalStructure.setStructureFormat(StructureFormatType.SMILES);
                    chemicalStructure.setStructureData(value);
                }
                reference.setStructure(chemicalStructure);
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

        InternationalBrandListType brands = drug.getInternationalBrands();
        if(brands != null) {
            for (InternationalBrandType brandType : brands.getInternationalBrand()) {
                named.addName(brandType.getName());
            }
        }
    }

    private String createRDFId(DrugType drug) {
        String primaryDrugId = getPrimaryDrugId(drug);
        return primaryDrugId != null ? primaryDrugId : UUID.randomUUID() + "";
    }

    private String getPrimaryDrugId(DrugType drug) {
        for (DrugbankDrugSaltIdType drugbankDrugIdType : drug.getDrugbankId()) {
            if(drugbankDrugIdType.isPrimary()) {
                return drugbankDrugIdType.getValue();
            }
        }

        return null;
    }

}
