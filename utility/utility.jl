mutable struct StopWatch
    start::Float
    interval::Float
    f::Function
    StopWatch(_interval, callback) = new(time(), _interval, callback)
end

function check(watch::StopWatch)
    now = time()
    if now - watch.start > watch.interval
        watch.f()
        watch.start = now
    end
end

function progressBar(step, total)
    barWidth = 70
    percent = round(step / total * 100.0, digits = 2)
    str = "["
    pos = barWidth * percent / 100.0
    for i in 1:barWidth
        if i <= pos
            str *= "I"
        else
            str *= " "
        end
    end
    str *= "] $step/$total=$percent%"
    return str
end
