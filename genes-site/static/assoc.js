'use strict';

LocusZoom.Data.GwasCatalog.prototype.combineChainBody_orig = LocusZoom.Data.GwasCatalog.prototype.combineChainBody;
LocusZoom.Data.GwasCatalog.prototype.combineChainBody = function (data, chain, fields, outnames, trans) {
    // this detects if we are preparing data for the annotation_catalog track
    // if we are, then it returns all of the GWASCatalog data for this region.
    // otherwise, it does the normal behavior of returning the association variants leftmerged with the GWASCatalog data.
    // (ie, normally GWASCatalog hits are only shown if they overlap the position of one of the variants being plotted)
    if (data.length && fields.length == 4 && fields.indexOf('variant')!==-1) {
        //console.log({data:data, chain:chain, fields:fields, outnames:outnames, trans:trans});
        return data.map(function(d) {
            return {
                'assoc:variant': d.variant, 'assoc:chromosome': d.chrom, 'assoc:position': d.pos,
                'catalog:variant': d.variant, 'catalog:rsid': d.rsid, 'catalog:trait': d.trait, 'catalog:log_pvalue': d.log_pvalue
            }
        });
    }
    return LocusZoom.Data.GwasCatalog.prototype.combineChainBody_orig(data, chain, fields, outnames, trans);
};

$.getJSON('/api/variants/'+model.genename+'/'+model.phecode).then(function(resp) {
    window._debug.resp = resp;
    var chrom = resp.chrom;

    var df = resp.df;
    df.pos = df.pos_delta;
    delete df.pos_delta;
    _.map(df.pos, function(elem, i) { if (i>0) { df.pos[i] += df.pos[i-1]; } })
    var data = _.sortBy(dataframe_to_objects(df), _.property('pos'));

    // Plot
    var variants = objects_to_dataframe(data);
    variants.chromosome = variants.pos.map(function() {return chrom;});
    variants.id = variants.pos.map(function(pos, i) { return fmt('{0}:{1}_{2}', chrom, pos, variants.base[i].replace(':','/')); });
    variants.variant = variants.id;
    variants.pvalue = variants.pval; delete variants.pval;
    variants.position = variants.pos; delete variants.pos;
    //variants.log_pvalue = variants.pval.map(function(p) { return -Math.log10(Math.max(1e-100, p)); });

    var start_position = d3.min(variants.position);
    var end_position = d3.max(variants.position);

    LocusZoom.TransformationFunctions.set("neglog10_or_100", function(x){return (x === 0) ? 100 : -Math.log10(x);});

    var remoteBase = "https://portaldev.sph.umich.edu/api/v1/";
    var data_sources = new LocusZoom.DataSources()
        .add("assoc", ["StaticJSON", variants])
        .add("catalog", ["GwasCatalogLZ", {url: remoteBase + 'annotation/gwascatalog/results/'}])
        .add("ld", ["LDLZ", {url: remoteBase + "pair/LD/", params: { pvalue_field: "assoc:pvalue|neglog10_or_100" }}])
        .add("gene", ["GeneLZ", { url: remoteBase + "annotation/genes/"}])
        .add("recomb", ["RecombLZ", { url: remoteBase + "annotation/recomb/results/" }]);

    var layout = {
        unnamedspaced: true,
        width: 800,
        height: 400,
        max_region_scale: 1e7, // arbitrary but high enough to show all data
        responsive_resize: 'width_only',
        mouse_guide: false,
        dashboard: {components: [ {type: 'download', position: 'right', color: 'gray' } ]},
        state: {
            chr: chrom,
            start: Math.round(start_position - 10 - (end_position-start_position)/30),
            end: Math.round(end_position + 10 + (end_position-start_position)/30),
            genome_build: 'GRCh38',
        },
        panels: [
            custom_LocusZoom_Layouts_get("panel", "annotation_catalog", {
                unnamespaced: true,
                height: 52, min_height: 52,
                margin: { top: 30, bottom: 13 },
                dashboard: { components: [] },
                title: {
                    text: 'Hits in GWAS Catalog',
                    style: {'font-size': '14px'},
                    x: 50,
                },
                "data_layers.0.hit_area_width": 50,
            }),
            custom_LocusZoom_Layouts_get("panel", "association_catalog", {
                unnamedspaced: true,
                height: 200, min_height: 200,
                dashboard: { components: [{type:'toggle_legend', position:'right'}] },
                data_layers: [
                    LocusZoom.Layouts.get("data_layer", "recomb_rate", { unnamespaced: true }),
                    LocusZoom.Layouts.get("data_layer", "association_pvalues_catalog", {
                        unnamespaced: true,
                        fields: [
                            "{{namespace[assoc]}}id",
                            "{{namespace[assoc]}}position",
                            "{{namespace[assoc]}}pvalue|neglog10_or_100",
                            "{{namespace[assoc]}}pvalue|scinotation",
                            "{{namespace[assoc]}}maf|scinotation", "{{namespace[assoc]}}mac_case", "{{namespace[assoc]}}mac_control",
                            "{{namespace[ld]}}state", "{{namespace[ld]}}isrefvar",
                            "{{namespace[catalog]}}rsid", "{{namespace[catalog]}}trait", "{{namespace[catalog]}}log_pvalue"
                        ],
                        id_field: "{{namespace[assoc]}}id",
                        tooltip: {
                            closable: true,
                            show: {
                                "or": ["highlighted", "selected"]
                            },
                            hide: {
                                "and": ["unhighlighted", "unselected"]
                            },
                            html: "<strong>chr{{{{namespace[assoc]}}id}}</strong><br>" +
                                "<span style='white-space:nowrap'>P-value: <strong>{{{{namespace[assoc]}}pvalue|scinotation}}</strong></span><br>" +
                                "MAF: <strong>{{{{namespace[assoc]}}maf|scinotation}}</strong><br>" +
                                "Case MAC: <strong>{{{{namespace[assoc]}}mac_case}}</strong><br>" +
                                "Control MAC: <strong>{{{{namespace[assoc]}}mac_control}}</strong><br>" +
                                "{{#if {{namespace[catalog]}}rsid}}<br><a href=\"https://www.ebi.ac.uk/gwas/search?query={{{{namespace[catalog]}}rsid}}\" target=\"_new\">See hits in GWAS catalog</a>{{/if}}" +
                                "<br><a href=\"javascript:void(0);\" onclick=\"LocusZoom.getToolTipDataLayer(this).makeLDReference(LocusZoom.getToolTipData(this));\">Make LD Reference</a>"
                        },
                        x_axis: { field: "{{namespace[assoc]}}position" },
                        y_axis: {
                            axis: 1,
                            field: "{{namespace[assoc]}}pvalue|neglog10_or_100",
                            floor: 0,
                            upper_buffer: 0.1,
                            min_extent: [0, 10]
                        }
                    })
                ],
                "legend.origin.y": 15,
                "margin.top": 10,
            }),
            custom_LocusZoom_Layouts_get("panel", "genes", {
                unnamespaced: true,
                dashboard: {components: []},
                data_layers: [
                    LocusZoom.Layouts.get("data_layer", "genes", {
                        unnamespaced: true,
                        fields: ["{{namespace[gene]}}all"],
                        tooltip: {
                            closable: true,
                            show: {
                                or: ["highlighted", "selected"]
                            },
                            hide: {
                                and: ["unhighlighted", "unselected"]
                            },
                            html: "<h4><strong><i>{{gene_name}}</i></strong></h4><div>Gene ID: <strong>{{gene_id}}</strong></div><div>Transcript ID: <strong>{{transcript_id}}</strong></div><div style=\"clear: both;\"></div><table width=\"100%\"><tr><td style=\"text-align: right;\"><a href=\"http://exac.broadinstitute.org/gene/{{gene_id}}\" target=\"_new\">More data on ExAC</a></td></tr></table>"
                        },
                        label_exon_spacing: 3,
                        exon_height: 8,
                        bounding_box_padding: 5,
                        track_vertical_spacing: 5
                    })
                ],
                "margin.top": 0,
            })
        ]
    };
    LocusZoom.Layouts.add("plot", "custom_association", layout);
    layout = LocusZoom.Layouts.get("plot", "custom_association");

    window._debug.variants = variants;
    $(function() {
        var plot = LocusZoom.populate("#lz_plot_container", data_sources, layout);
        // resize the genes panel to fit the data, but only once.
        // Resizing the genes panel is disabled for now due to bugs.  Sometimes the genes data_layer ends up inside the association panel.
        // TODO: give a max-height of 500px.
        // var gene_resize_hook = plot.panels.genes.on('data_rendered', function(){
        //     plot.panels.genes.scaleHeightToData();
        //     plot.panels.genes.off('data_rendered', gene_resize_hook);
        // });
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
                {title: 'Position on chr'+resp.chrom, field:'pos', formatter:'comma_fmt'},
                {title: 'Allele', field:'base'},
                {title: 'MAF (Minor Allele Frequency)', field:'maf', formatter:'2digit_fmt'},
                {title: 'Case MAC (Minor Allele Count)', field:'mac_case'},
                {title: 'Control MAC (Minor Allele Count)', field:'mac_control'},
                {title: 'P-value', field:'pval', formatter:'2digit_fmt'},
            ],
            data: data,
            initialSort: [{column:'pval', dir:'asc'}],
            tooltipGenerationMode:'hover',tooltips:tabulator_tooltip_maker,tooltipsHeader:true,
        });
    });
});
