//
// This file was generated by the JavaTM Architecture for XML Binding(JAXB) Reference Implementation, v2.2.7 
// See <a href="http://java.sun.com/xml/jaxb">http://java.sun.com/xml/jaxb</a> 
// Any modifications to this file will be lost upon recompilation of the source schema. 
// Generated on: 2017.02.06 at 04:04:57 PM EST 
//


package ca.drugbank.model;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.List;
import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlType;


/**
 * <p>Java class for carrier-type complex type.
 * 
 * <p>The following schema fragment specifies the expected content contained within this class.
 * 
 * <pre>
 * &lt;complexType name="carrier-type">
 *   &lt;complexContent>
 *     &lt;restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
 *       &lt;sequence>
 *         &lt;group ref="{http://www.drugbank.ca}interactant-group"/>
 *       &lt;/sequence>
 *       &lt;attribute name="position" type="{http://www.w3.org/2001/XMLSchema}integer" />
 *     &lt;/restriction>
 *   &lt;/complexContent>
 * &lt;/complexType>
 * </pre>
 * 
 * 
 */
@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "carrier-type", propOrder = {
    "id",
    "name",
    "organism",
    "actions",
    "references",
    "knownAction",
    "polypeptide"
})
public class CarrierType {

    @XmlElement(required = true)
    protected String id;
    @XmlElement(required = true)
    protected String name;
    @XmlElement(required = true)
    protected String organism;
    @XmlElement(required = true)
    protected ActionListType actions;
    @XmlElement(required = true)
    protected ReferenceListType references;
    @XmlElement(name = "known-action", required = true)
    protected KnownActionType knownAction;
    protected List<PolypeptideType> polypeptide;
    @XmlAttribute(name = "position")
    protected BigInteger position;

    /**
     * Gets the value of the id property.
     * 
     * @return
     *     possible object is
     *     {@link String }
     *     
     */
    public String getId() {
        return id;
    }

    /**
     * Sets the value of the id property.
     * 
     * @param value
     *     allowed object is
     *     {@link String }
     *     
     */
    public void setId(String value) {
        this.id = value;
    }

    /**
     * Gets the value of the name property.
     * 
     * @return
     *     possible object is
     *     {@link String }
     *     
     */
    public String getName() {
        return name;
    }

    /**
     * Sets the value of the name property.
     * 
     * @param value
     *     allowed object is
     *     {@link String }
     *     
     */
    public void setName(String value) {
        this.name = value;
    }

    /**
     * Gets the value of the organism property.
     * 
     * @return
     *     possible object is
     *     {@link String }
     *     
     */
    public String getOrganism() {
        return organism;
    }

    /**
     * Sets the value of the organism property.
     * 
     * @param value
     *     allowed object is
     *     {@link String }
     *     
     */
    public void setOrganism(String value) {
        this.organism = value;
    }

    /**
     * Gets the value of the actions property.
     * 
     * @return
     *     possible object is
     *     {@link ActionListType }
     *     
     */
    public ActionListType getActions() {
        return actions;
    }

    /**
     * Sets the value of the actions property.
     * 
     * @param value
     *     allowed object is
     *     {@link ActionListType }
     *     
     */
    public void setActions(ActionListType value) {
        this.actions = value;
    }

    /**
     * Gets the value of the references property.
     * 
     * @return
     *     possible object is
     *     {@link ReferenceListType }
     *     
     */
    public ReferenceListType getReferences() {
        return references;
    }

    /**
     * Sets the value of the references property.
     * 
     * @param value
     *     allowed object is
     *     {@link ReferenceListType }
     *     
     */
    public void setReferences(ReferenceListType value) {
        this.references = value;
    }

    /**
     * Gets the value of the knownAction property.
     * 
     * @return
     *     possible object is
     *     {@link KnownActionType }
     *     
     */
    public KnownActionType getKnownAction() {
        return knownAction;
    }

    /**
     * Sets the value of the knownAction property.
     * 
     * @param value
     *     allowed object is
     *     {@link KnownActionType }
     *     
     */
    public void setKnownAction(KnownActionType value) {
        this.knownAction = value;
    }

    /**
     * Gets the value of the polypeptide property.
     * 
     * <p>
     * This accessor method returns a reference to the live list,
     * not a snapshot. Therefore any modification you make to the
     * returned list will be present inside the JAXB object.
     * This is why there is not a <CODE>set</CODE> method for the polypeptide property.
     * 
     * <p>
     * For example, to add a new item, do as follows:
     * <pre>
     *    getPolypeptide().add(newItem);
     * </pre>
     * 
     * 
     * <p>
     * Objects of the following type(s) are allowed in the list
     * {@link PolypeptideType }
     * 
     * 
     */
    public List<PolypeptideType> getPolypeptide() {
        if (polypeptide == null) {
            polypeptide = new ArrayList<PolypeptideType>();
        }
        return this.polypeptide;
    }

    /**
     * Gets the value of the position property.
     * 
     * @return
     *     possible object is
     *     {@link BigInteger }
     *     
     */
    public BigInteger getPosition() {
        return position;
    }

    /**
     * Sets the value of the position property.
     * 
     * @param value
     *     allowed object is
     *     {@link BigInteger }
     *     
     */
    public void setPosition(BigInteger value) {
        this.position = value;
    }

}
