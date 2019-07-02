'use strict';

LocusZoom.TransformationFunctions.add('commas', function(x){return x.toLocaleString()});

$.getJSON('/api/gene/'+model.genename).then(function(resp) {
    // fields = phecode, phenostring, category, pval
    var assocs = resp.assocs;
    assocs.id = assocs.phecode;
    assocs.trait_label = assocs.phecode.map(function(d, i) { return assocs.phecode[i] + ' - ' + assocs.phenostring[i]; });;
    assocs.log_pvalue = assocs.pval.map(function(p) { return -Math.log10(Math.max(1e-100, p)); });
    assocs.trait_group = assocs.category;

    var significance_threshold = -Math.log10(0.05 / assocs.id.length);
    var best_nlpval = d3.max(assocs.log_pvalue);

    var data_sources = new LocusZoom.DataSources().add('phewas', ['StaticJSON', assocs]);
    var layout = {
        width: 800,
        height: 400,
        min_width: 800,
        min_height: 400,
        responsive_resize: 'width_only',
        mouse_guide: false,
        dashboard: {components: [ {type: 'download', position: 'right', color: 'gray' } ]},
        panels: [
            LocusZoom.Layouts.get('panel', 'phewas', {
                margin: {top: 5, right: 5, bottom: 120, left: 50 }
            })
        ],
    }

    layout.panels[0].data_layers[0].offset = significance_threshold;
    layout.panels[0].data_layers[1].fields.push('phewas:phenostring');
    layout.panels[0].data_layers[1].fields.push('phewas:num_rare');
    layout.panels[0].data_layers[1].fields.push('phewas:startpos', 'phewas:endpos');
    layout.panels[0].data_layers[1].fields.push('phewas:mac_case', 'phewas:mac_control');
    layout.panels[0].data_layers[1].fields.push('phewas:num_cases', 'phewas:num_controls');
    layout.panels[0].data_layers[1].tooltip.html =
        ("<strong>{{phewas:trait_label|htmlescape}}</strong><br>" +
         "Category: <strong>{{phewas:trait_group|htmlescape}}</strong><br>" +
         "P-value: <strong>{{phewas:log_pvalue|logtoscinotation|htmlescape}}</strong><br>" +
         "#Rare Variants: <strong>{{phewas:num_rare}}</strong><br>" +
         "#Cases / #Controls: <strong>{{phewas:num_cases|commas}} / {{phewas:num_controls|commas}}</strong><br>" +
         "Case / Control MAC: <strong>{{phewas:mac_case}} / {{phewas:mac_control}}</strong><br>" +
         "Start - End: <strong>{{phewas:startpos|commas}} - {{phewas:endpos|commas}}</strong><br>"
        );
    layout.panels[0].data_layers[1].behaviors.onclick = [{action: 'link', href: '/assoc/'+model.genename+'/{{phewas:id}}'}];
    layout.panels[0].data_layers[1].y_axis.min_extent = [0, significance_threshold*1.1];

    if (assocs.id.length <= 10) {
        layout.panels[0].data_layers[1].label.filters = []; // show all labels
    } else if (assocs.log_pvalue.filter(function(nlpval) { return nlpval == best_nlpval; }).length >= 6) {
        layout.panels[0].data_layers[1].label = false; // too many are tied for 1st and will make a mess so just hide all labels
    } else {
        var eighth_best_nlpval = _.sortBy(assocs.log_pvalue).reverse()[8];
        layout.panels[0].data_layers[1].label.filters = [
            {field: 'phewas:log_pvalue', operator: '>', value: significance_threshold},
            {field: 'phewas:log_pvalue', operator: '>', value: best_nlpval*0.5}, // must be in top half of screen
            {field: 'phewas:log_pvalue', operator: '>=', value: eighth_best_nlpval} // must be among the best
        ];
    }

    window._debug.assocs = assocs;
    $(function() {
        var plot = LocusZoom.populate("#phewas_plot_container", data_sources, layout);
        window._debug.plot = plot;
    });

    $(function() {
        var data = dataframe_to_objects(assocs);
        var table = new Tabulator('#table', {
            //height: 600, // setting height lets Tabulator's VirtualDOM load really fast but makes scrolling awkward
            layout: 'fitColumns',
            pagination: 'local',
            paginationSize: 15,
            columns: [
                {title: 'Category', field:'category', headerFilter:true, widthGrow:1},
                {title: 'Code', field:'phecode', formatter:'link', formatterParams: { urlPrefix: '/assoc/'+model.genename+'/' }, headerFilter:true},
                {title: 'Name', field:'phenostring', formatter:'link', formatterParams: { url: function(cell){return '/assoc/'+model.genename+'/'+cell.getData().phecode;}}, headerFilter:true, widthGrow:2},
                {title: '#Cases', field:'num_cases', formatter:'comma_fmt'},
                {title: '#Controls', field:'num_controls', formatter:'comma_fmt'},
                {title: 'P-value', field:'pval', formatter:'2digit_fmt'},
                {title: '#Rare Variants', field:'num_rare', formatter:'comma_fmt'},
                {title: 'Start-End', field:'startpos', formatter: function(cell){return cell.getValue().toLocaleString()+' - '+cell.getData().endpos.toLocaleString()}, widthGrow:1},
                {title: 'Case MAC (Minor Allele Count)', field:'mac_case'},
                {title: 'Control MAC (Minor Allele Count)', field:'mac_control'},
            ],
            data: data,
            initialSort: [{column:'pval', dir:'asc'}],
            tooltipGenerationMode:'hover',tooltips:tabulator_tooltip_maker,tooltipsHeader:true,
        });
    });
});
