'use strict';

$.getJSON('/api/variants/'+model.genename+'/'+model.phecode).then(function(resp) {
    window._debug.resp = resp;
    var chrom = resp.chrom;

    var df = resp.df;
    df.pos = df.pos_delta;
    delete df.pos_delta;
    _.map(df.pos, function(elem, i) { if (i>0) { df.pos[i] += df.pos[i-1]; } })
    var data = _.sortBy(dataframe_to_objects(df), _.property('pos'));

    // Plot
    // TODO: switch to a normal LZ plot
    var variants = objects_to_dataframe(data);
    variants.id = variants.pos.map(function(pos, i) {
        return fmt('{0}-{1}-{2}', chrom, pos, variants.base[i]);
    });
    variants.trait_label = variants.id;
    variants.trait_group = variants.pos.map(function() { return 'chr'+chrom });
    variants.log_pvalue = variants.pval.map(function(p) { return -Math.log10(Math.max(1e-100, p)); });

    var significance_threshold = 3;
    var best_nlpval = d3.max(variants.log_pvalue);

    var data_sources = new LocusZoom.DataSources().add('phewas', ['StaticJSON', variants]);
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
                margin: {top: 5, right: 5, bottom: 50, left: 50 }
            })
        ],
    };

    layout.panels[0].data_layers[0].offset = significance_threshold;
    layout.panels[0].data_layers[1].fields.push('phewas:pos');
    layout.panels[0].data_layers[1].fields.push('phewas:maf');
    layout.panels[0].data_layers[1].fields.push('phewas:mac_case', 'phewas:mac_control');
    layout.panels[0].data_layers[1].tooltip.html =
        ("<strong>chr{{phewas:trait_label|htmlescape}}</strong><br>" +
         "P-value: <strong>{{phewas:log_pvalue|logtoscinotation|htmlescape}}</strong><br>" +
         "Case / Control MAC: <strong>{{phewas:mac_case}} / {{phewas:mac_control}}</strong><br>" +
         "MAF: <strong>{{phewas:maf}}</strong>"
        );
    layout.panels[0].data_layers[1].y_axis.min_extent = [0, significance_threshold*1.1];

    if (variants.id.length <= 10) {
        layout.panels[0].data_layers[1].label.filters = []; // show all labels
    } else if (variants.log_pvalue.filter(function(nlpval) { return nlpval == best_nlpval; }).length >= 6) {
        layout.panels[0].data_layers[1].label = false; // too many are tied for 1st and will make a mess so just hide all labels
    } else {
        var eighth_best_nlpval = _.sortBy(variants.log_pvalue).reverse()[8];
        layout.panels[0].data_layers[1].label.filters = [
            {field: 'phewas:log_pvalue', operator: '>', value: best_nlpval*0.5}, // must be in top half of screen
            {field: 'phewas:log_pvalue', operator: '>=', value: eighth_best_nlpval} // must be among the best
        ];
    }

    window._debug.variants = variants;
    $(function() {
        var plot = LocusZoom.populate("#phewas_plot_container", data_sources, layout);
        window._debug.plot = plot;
    });


    // Table
    $(function() {
        var table = new Tabulator('#table', {
            //height: 600, // setting height lets Tabulator's VirtualDOM load really fast but makes scrolling awkward
            layout: 'fitColumns',
            pagination: 'local', // TODO: disable pagination if <100 variants
            paginationSize: 100,
            columns: [
                {title: 'Position on chr'+resp.chrom, field:'pos'},
                {title: 'Allele', field:'base'},
                {title: 'MAF (Minor Allele Frequency)', field:'maf'},
                {title: 'Case MAC (Minor Allele Count)', field:'mac_case'},
                {title: 'Control MAC (Minor Allele Count)', field:'mac_control'},
                {title: 'P-value', field:'pval'},
            ],
            data: data,
            initialSort: [{column:'pval', dir:'asc'}],
            tooltipGenerationMode:'hover',tooltips:tabulator_tooltip_maker,tooltipsHeader:true,
        });
    });
});
