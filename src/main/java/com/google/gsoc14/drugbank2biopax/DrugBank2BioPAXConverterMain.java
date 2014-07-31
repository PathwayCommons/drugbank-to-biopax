package com.google.gsoc14.drugbank2biopax;

import com.google.gsoc14.drugbank2biopax.controller.DrugBank2BioPAXConverter;
import org.apache.commons.cli.*;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.biopax.paxtools.io.SimpleIOHandler;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.trove.TProvider;
import org.biopax.paxtools.util.BPCollections;

import javax.xml.bind.JAXBException;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;

public class DrugBank2BioPAXConverterMain {
    private static Log log = LogFactory.getLog(DrugBank2BioPAXConverterMain.class);
    private static final String helpText = DrugBank2BioPAXConverterMain.class.getSimpleName();

    public static void main( String[] args ) throws JAXBException {
        final CommandLineParser clParser = new GnuParser();
        Options gnuOptions = new Options();
        gnuOptions
                .addOption("d", "drugs", true, "structured DrugBank data (XML) [required]")
                .addOption("o", "output", true, "Output (BioPAX file) [required]")
        ;

        try {
            CommandLine commandLine = clParser.parse(gnuOptions, args);

            // DrugBank file and output file name are required!
            if(!commandLine.hasOption("o") && !commandLine.hasOption("d")) {
                HelpFormatter helpFormatter = new HelpFormatter();
                helpFormatter.printHelp(helpText, gnuOptions);
                System.exit(-1);
            }

            // Memory efficiency fix for huge BioPAX models
            BPCollections.I.setProvider(new TProvider());

            String drugBankFile = commandLine.getOptionValue("d");
            log.debug("Using DrugBank file: " + drugBankFile);
            FileInputStream drugBankStream = new FileInputStream(drugBankFile);

            DrugBank2BioPAXConverter drugBank2BioPAXConverter = new DrugBank2BioPAXConverter();
            Model model = drugBank2BioPAXConverter.convert(drugBankStream);

            String outputFile = commandLine.getOptionValue("o");
            log.debug("Saving the final BioPAX model to:" + outputFile);
            SimpleIOHandler simpleIOHandler = new SimpleIOHandler();
            FileOutputStream outputStream = new FileOutputStream(outputFile);
            simpleIOHandler.convertToOWL(model, outputStream);
            outputStream.close();
            log.debug("All done.");
        } catch (ParseException e) {
            System.err.println(e.getMessage());
            HelpFormatter helpFormatter = new HelpFormatter();
            helpFormatter.printHelp(helpText, gnuOptions);
            System.exit(-1);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }
}
