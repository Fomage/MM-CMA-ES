var cmaesFunction = function () { console.log("Warning: no function installed") }

var cmaesDefined

var cmaesXMean

var cmaesEigenVectors

var cmaesAxisLengths

var cmaesSteps

var cmaesPoints

var cmaesTimeout

var cmaesRedraw

var cmaesUpdateSlider = function() {}

var cmaesTerminationMessageElt

var colours = [["blue","violet"],["red","orange"],["green","teal"],["yellow","brown"]];

function cmaesPushStep() {
    cmaesSteps.push({
		defined: cmaesDefined,
        points: cmaesPoints,
        xmean: cmaesXMean,
        eigenVectors: cmaesEigenVectors,
        axisLengths: cmaesAxisLengths,
		sigma: cmaesSigma,
		//shouldSplit : cmaesShouldSplit,
		//divisionThreshold : cmaesDivisionThreshold
    })
    cmaesUpdateSlider()
}

function cmaesInitialize() {
    if (cmaesTimeout) {
        clearTimeout(cmaesTimeout)
        cmaesTimeout = null
    }
    _cmaesInitialize()
    cmaesLoop()
}

function cmaesLoop() {
    if (_cmaesStep()) {
        cmaesTimeout = setTimeout(cmaesLoop, 0)
    }
}

function cmaesSetTerminationMessage(msg) {
    while (cmaesTerminationMessageElt.hasChildNodes()) {
        cmaesTerminationMessageElt.removeChild(cmaesTerminationMessageElt.firstChild)
    }
    var text = document.createTextNode(msg)
    cmaesTerminationMessageElt.appendChild(text)
}

window.onload = function() {
    var canvasBox = document.getElementById("canvasBox")
    var canvas = document.getElementById("canvas")
    var formula = document.getElementById("formula")
    var form = document.getElementById("form")
    var drawButton = document.getElementById("drawButton")
    var slider = document.getElementById("slider")
    var terminationMessage = document.getElementById("terminationMessage")
    cmaesTerminationMessageElt = terminationMessage

    slider.valueAsNumber = 0

    updateFunction()
    function updateFunction() {
        eval("cmaesFunction = function (x, y) { return " + formula.value + " }")
        cmaesSteps = []
        cmaesInitialize()
    }

    var size
    window.onresize = function() {
        var width = canvasBox.clientWidth
        var height = canvasBox.clientHeight
        var newSize = Math.min(width, height)
        canvas.style.left = Math.floor((width - newSize) / 2) + "px"
        canvas.style.top = Math.floor((height - newSize) / 2) + "px"
        if (size == newSize) {
            return
        }
        size = newSize
        canvas.width = size
        canvas.height = size
        redraw()
    }
    window.onresize()

    form.onsubmit = function(ev) {
        updateFunction()
        redraw()
        ev.preventDefault()
    }

    function canvasOfPlan(planPoint) {
        return [
            Math.round((planPoint[0] + 0.5) * size),
            Math.round((planPoint[1] + 0.5) * size)]
    }

    function planOfCanvas(canvasPoint) {
        return [
            canvasPoint[0] / size - 0.5,
            canvasPoint[1] / size - 0.5]
    }

    function translatePoint(point, translation) {
        return [
            point[0] + translation[0],
            point[1] + translation[1]]
    }

    function redraw() {
        var ctx = canvas.getContext("2d")    
        var imageData = ctx.createImageData(size, size)
        var data = imageData.data
        for (var y = 0; y < size; y++) {
            for (var x = 0; x < size; x++) {
                var plan = planOfCanvas([x, y])
                var v = Math.round(cmaesFunction(plan[0], plan[1]) * 256)
                var pixelIndex = (y * size + x) * 4
                data[pixelIndex] = v
                data[pixelIndex + 1] = v
                data[pixelIndex + 2] = v
                data[pixelIndex + 3] = 255
            }
        }
        ctx.putImageData(imageData, 0, 0)
        var step = cmaesSteps[slider.valueAsNumber]
        if (step) {
			var pointRank=0;
			for(ivillage=0;ivillage<step.xmean.length;ivillage++){if(step.defined[ivillage]){
				step.points[ivillage].forEach(function (planPoint) {//draws distribution
					var canvasPoint = canvasOfPlan([planPoint[0],planPoint[1]])
					var radius = 2
					/*if(planPoint[2]!=0){
						radius=4
					}*/
					ctx.beginPath()
					ctx.arc(canvasPoint[0], canvasPoint[1], radius, 0, 2 * Math.PI, false)
					/*if(step.shouldSplit!=0){
						if(planPoint[3]==1){
							ctx.fillStyle = "green"
							ctx.strokeStyle = "blue"
						}else{
							ctx.fillStyle = "yellow"
							ctx.strokeStyle = "pink"
						}
					}else{
						ctx.fillStyle = "black"
						if(pointRank<2){
							ctx.strokeStyle = "red"
						}else{
							ctx.strokeStyle = "white"
						}
					}*/
					ctx.fillStyle=colours[ivillage][0]
					ctx.strokeStyle=colours[ivillage][1]
					ctx.fill()
					ctx.stroke()
					pointRank++
				})
				var canvasXmean = canvasOfPlan(step.xmean[ivillage])//draws xmean
				var path = new Path2D()
				var radius = 4
				ctx.beginPath()
				ctx.arc(canvasXmean[0], canvasXmean[1], radius, 0, 2 * Math.PI, false)
				/*if(step.shouldSplit[ivillage]!=0){
					ctx.fillStyle = "red"
				}else{
					ctx.fillStyle = "yellow"
				}*/
				ctx.fillStyle = colours[ivillage][0]
				ctx.strokeStyle = colours[ivillage][1]
				ctx.fill()
				ctx.stroke()
				
				/*var scaleFacEigenvactors = 1/4//draws eigenvector
				//draws first point
				var eigenV1 = canvasOfPlan([step.eigenVectors[0][0]*scaleFacEigenvactors+step.xmean[0],
					step.eigenVectors[0][1]*scaleFacEigenvactors+step.xmean[1]])
				path = new Path2D()
				radius = 4
				ctx.beginPath()
				ctx.arc(eigenV1[0], eigenV1[1], radius, 0, 2 * Math.PI, false)
				ctx.fillStyle = "blue"
				ctx.strokeStyle = "orange"
				ctx.fill()
				ctx.stroke()
				//draws second point
				var eigenV2 = canvasOfPlan([step.eigenVectors[1][0]*scaleFacEigenvactors+step.xmean[0],
					step.eigenVectors[1][1]*scaleFacEigenvactors+step.xmean[1]])
				path = new Path2D()
				radius = 4
				ctx.beginPath()
				ctx.arc(eigenV2[0], eigenV2[1], radius, 0, 2 * Math.PI, false)
				ctx.fillStyle = "blue"
				ctx.strokeStyle = "orange"
				ctx.fill()
				ctx.stroke()
				//draws lines
				ctx.beginPath()
				ctx.moveTo(eigenV1[0],eigenV1[1])
				ctx.lineTo(canvasXmean[0],canvasXmean[1])
				ctx.lineTo(eigenV2[0],eigenV2[1])
				ctx.strokeStyle = "blue"
				ctx.stroke()*/
				
				/*ctx.beginPath()//draws division threshold
				var x = canvasXmean[0]
				var y = canvasXmean[1]
				var divisionBoundaryRadius=step.divisionThreshold
				ctx.arc(x,y,divisionBoundaryRadius*3*size/2,0,2*Math.PI)
				ctx.strokeStyle = "purple"
				ctx.stroke()*/

				ctx.beginPath()//draws ellipse
				var scaleFacEllipse = 2
				var rx = step.axisLengths[ivillage][1] * size * step.sigma[ivillage] * scaleFacEllipse / 2
				var ry = step.axisLengths[ivillage][0] * size * step.sigma[ivillage] * scaleFacEllipse / 2
				var angle = Math.atan2(
					step.eigenVectors[ivillage][0][0],
					step.eigenVectors[ivillage][0][1])// * 180 / Math.PI
				/*cmaesSetTerminationMessage("Sigma="+step.sigma[ivillage]+" xmean="+step.xmean[ivillage]+" f(xmean)="
					+cmaesFunction(step.xmean[ivillage][0],step.xmean[ivillage][1])
					+" divisionThreshold="+step.divisionThreshold[ivillage])*/
				ctx.ellipse(canvasXmean[0], canvasXmean[1], rx, ry, -angle, 0, 2 * Math.PI, false)
				ctx.strokeStyle = colours[ivillage][1]
				ctx.stroke()
			}}
            
        }
    }
    cmaesRedraw = redraw

    function updateSlider() {
        var position = slider.valueAsNumber
        var onmax = position == slider.max
        slider.max = cmaesSteps.length - 1
        if (onmax) {
            slider.valueAsNumber = slider.max
            redraw()
        }
    }
    cmaesUpdateSlider = updateSlider

    slider.onchange = slider.oninput = function() {
        redraw()
    }
}

