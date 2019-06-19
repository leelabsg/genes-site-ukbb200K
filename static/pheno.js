'use strict';

$.getJSON('/api/pheno/'+model.phecode).then(function(resp) {
    var num_genes = resp.num_genes;
    var assocs = objects_to_dataframe(_.sortBy(dataframe_to_objects(resp.assocs), _.property('startpos')));

    assocs.id = assocs.name;
    assocs.trait_label = assocs.name;
    assocs.log_pvalue = assocs.pval.map(function(p) { return -Math.log10(Math.max(1e-6, p)); });
    assocs.trait_group = assocs.chrom.map(function(s) { return 'chr'+s.toString().padStart(2,'0'); });

    var significance_threshold = -Math.log10(0.05 / num_genes);
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
                margin: {top: 5, right: 5, bottom: 110, left: 50 }
            })
        ],
    }

    layout.panels[0].data_layers[0].offset = significance_threshold;
    layout.panels[0].data_layers[1].fields.push('phewas:num_rare');
    layout.panels[0].data_layers[1].fields.push('phewas:chrom', 'phewas:startpos', 'phewas:endpos');
    layout.panels[0].data_layers[1].fields.push('phewas:mac_case', 'phewas:mac_control');
    layout.panels[0].data_layers[1].tooltip.html =
        ("gene: <strong>{{phewas:trait_label|htmlescape}}</strong><br>" +
         "P-value: <strong>{{phewas:log_pvalue|logtoscinotation|htmlescape}}</strong><br>" +
         "#Rare Variants: <strong>{{phewas:num_rare}}</strong><br>" +
         "Case / Control MAC: <strong>{{phewas:mac_case}} / {{phewas:mac_control}}</strong><br>" +
         "Chrom:Start-End: <strong>{{phewas:chrom}} : {{phewas:startpos}} - {{phewas:endpos}}</strong><br>"
        );
    layout.panels[0].data_layers[1].behaviors.onclick = [{action: 'link', href: '/assoc/{{phewas:id}}/'+model.phecode}];
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
                {title: 'Gene', field:'name', formatter:'link', formatterParams: {url:function(cell){return '/assoc/'+cell.getValue()+'/'+model.phecode}}, headerFilter:true, widthGrow:3},
                {title: 'P-value', field:'pval'},
                {title: '#Rare Variants', field:'num_rare'},
                {title: 'Chromosome', field:'chrom', headerFilter:true, headerFilterFunc:'='},
                {title: 'Start-End', field:'startpos', formatter: function(cell){return cell.getValue()+' - '+cell.getData().endpos}},
                {title: 'Case MAC (Minor Allele Count)', field:'mac_case'},
                {title: 'Control MAC (Minor Allele Count)', field:'mac_control'},
            ],
            data: data,
            initialSort: [{column:'pval', dir:'asc'}],
            tooltipGenerationMode: 'hover', // generate tooltips just-in-time when the data is hovered
            tooltips: function(cell) {
                // this function attempts to check whether an ellipsis ('...') is hiding part of the data.
                // to do so, I compare element.clientWidth against element.scrollWidth;
                // when scrollWidth is bigger, that means we're hiding part of the data.
                // unfortunately, the ellipsis sometimes activates too early, meaning that even with clientWidth == scrollWidth some data is hidden by the ellipsis.
                // fortunately, these tooltips are just a convenience so I don't mind if they fail to show.
                // I don't know whether clientWidth or offsetWidth is better. clientWidth was more convenient in Chrome74.
                var e = cell.getElement();
                //return '' + e.offsetWidth + ' || ' + e.scrollWidth + ' || ' + e.clientWidth;
                if (e.clientWidth >= e.scrollWidth) {
                    return false; // all the text is shown, so there is no '...', so no tooltip is needed
                } else {
                    return cell.getValue();
                }
            }
        });
    });
});
