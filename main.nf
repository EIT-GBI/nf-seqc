// TODO: include processes from submodules under modules/
// include { TOOL } from './modules/<tool>/main.nf'

workflow {
    // TODO: build the input channel, e.g.:
    // ch_input = Channel.fromPath(params.input)
    //     .map { f -> tuple([id: f.baseName], f) }

    // TODO: wire processes, e.g.:
    // ch_tool = TOOL(ch_input)

    // Once you have channels to publish, split this into `main:` / `publish:`:
    //
    //     main:
    //     ch_tool = TOOL(ch_input)
    //
    //     publish:
    //     tool = ch_tool
    //
    // and uncomment the matching `output {}` block below.
}

// output {
//     tool {
//         path 'tool'
//     }
// }
