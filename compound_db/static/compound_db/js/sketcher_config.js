/**
 * Created by sichom on 3/18/16.
 */

var initChemSketcher = function() {
    var sketcher_height = 400;
    var options = {
        'oneMolecule' : true
        , 'useServices' : false
        , 'includeQuery' : true
    };
    var chem_sketcher = new ChemDoodle.SketcherCanvas(
        "chemdoodle_sketcher"
        , 500
        , sketcher_height
        , options
    );

    // sets terminal carbon labels to display
    chem_sketcher.specs.atoms_displayTerminalCarbonLabels_2D = true;
    // sets atom labels to be colored by JMol colors, which are easy to recognize
    chem_sketcher.specs.atoms_useJMOLColors = true;
    // enables overlap clear widths, so that some depth is introduced to overlapping bonds
    chem_sketcher.specs.bonds_clearOverlaps_2D = true;
    // sets the shape color to improve contrast when drawing figures
    chem_sketcher.specs.shapes_color = 'c10000';

    // resize the canvas according to wrapper size on window size change
    var resizeCanvas = function() {
        var sketcher = $(".sketcher");
        chem_sketcher.resize(sketcher.width(), sketcher_height)
    };
    $(window).resize(resizeCanvas);

    //
    var submitMol = function() {
        var mol = chem_sketcher.getMolecule();
        $(".mol_file_field").val(ChemDoodle.writeMOL(mol))
    };
    $(".submit-form-with-mol").click(submitMol);

    // make
    resizeCanvas();

    // force repaint - just to be sure
    chem_sketcher.repaint();
};